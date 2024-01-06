#include <filesystem>
#include <fstream>

#include <boost/program_options.hpp>
#include <omp.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spoa/spoa.hpp>

#include "thesis/corrector/fragmented_corrector.hpp"
#include "thesis/format/paf.hpp"
#include "thesis/utility/file_io/parse.hpp"
#include "thesis/utility/file_io/read_record.hpp"

namespace bpo = boost::program_options;
namespace fs = std::filesystem;

/* --- print some error message --- */
template <class T>
concept printable = requires(std::ostream& os, const T& obj) {
  { os << obj } -> std::same_as<std::ostream&>;
};

template <printable T = std::string>
void exit_and_print_help(T msg = ""s) {
  std::cerr << msg << std::endl;
  exit(EXIT_FAILURE);
}

auto make_path_checker(const std::string& opt_name) {
  return [opt_name](const fs::path& p) {
    if (!fs::exists(p)) {
      throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                  opt_name, p);
    }
  };
};

auto check_argument(const bpo::variables_map& vmap) {
  auto reads_path = vmap["raw_read"].as<fs::path>();
  auto overlap_path = vmap["overlap_info"].as<fs::path>();
  auto output_path = vmap["output_path"].as<fs::path>();

  /* check file format */
  check_file_format("raw_read", reads_path,
                    FILE_FORMAT::FASTA | FILE_FORMAT::FASTQ);
  check_file_format("overlap_info", overlap_path, FILE_FORMAT::PAF);
  auto fout = std::ofstream(output_path, std::ios::app);
  if (!fout.is_open()) {
    spdlog::error("Cannot open output file: {}", output_path.string());
    exit(1);
  }
  /* check platform */
  auto platform = vmap["platform"].as<std::string>();
  if (platform != "PacBio" && platform != "ONT") {
    spdlog::error("Platform must be PacBio or ONT");
    exit(1);
  }
  auto thread = vmap["thread"].as<int>();
  if (thread < 1) {
    spdlog::error("Thread must be greater than 0");
    exit(1);
  }
}

int main(int argc, char* argv[]) {
  /**
   * 1. raw reads
   * 2. overlap info
   * 3. output path
   * 4. sequencing platform
   */
  
  auto start = std::chrono::steady_clock::now();
  
  bpo::options_description opts{};
  bpo::variables_map vmap;
  try {
    opts.add_options()("help,h", "Show help message")(
        "raw_read,r",
        bpo::value<fs::path>()->required()->notifier(
            make_path_checker("raw_read")),
        "path to TGS raw reads")(
        "overlap_info,c",
        bpo::value<fs::path>()->required()->notifier(
            make_path_checker("overlap_info")),
        "path to the overlap information of raw reads(.paf)")(
        "output_path,o", bpo::value<fs::path>()->required(),
        "the output path of corrected read")(
        "platform,p", bpo::value<std::string>(),
        "the sequencing platform of raw reads, must be \"ONT\" or \"PacBio\"")(
        "thread,t", bpo::value<int>()->default_value(1),
        "the maximum number of threads")(
        "debug", bpo::bool_switch()->default_value(false), "enable debug mode");
    bpo::store(bpo::parse_command_line(argc, argv, opts), vmap);
    bpo::notify(vmap);
    if (vmap.contains("help")) {
      exit_and_print_help(opts);
    }
    check_argument(vmap);
  } catch (const std::exception& ex) {
    spdlog::error("{}", ex.what());
    exit_and_print_help(opts);
  }

  auto reads_path = vmap["raw_read"].as<fs::path>();
  auto overlap_path = vmap["overlap_info"].as<fs::path>();
  auto output_path = vmap["output_path"].as<fs::path>();
  auto platform = vmap["platform"].as<std::string>();

  auto debug_mode = vmap["debug"].as<bool>();
  if (debug_mode) {
    auto logger = spdlog::stderr_color_mt("logger");
    spdlog::set_default_logger(logger);
    spdlog::set_level(spdlog::level::debug);
  } else {
    spdlog::set_level(spdlog::level::info);
  }
  auto thread = vmap["thread"].as<int>();
  omp_set_num_threads(thread);

  auto overlap = read_records<bio::PafRecord>(overlap_path);
  auto call_corrector = [&]<class R>() {
    auto raw_reads = read_records<R>(reads_path);
    auto correcter = FragmentedReadCorrector<R>(std::move(raw_reads),
                                                std::move(overlap), platform,
                                                thread, debug_mode);
    auto corrected_read = correcter.correct();
    return corrected_read;
  };

  auto corrected_read = std::vector<bio::FastaRecord<false>>{};
  if (parse_file_format(reads_path) == FILE_FORMAT::FASTA) {
    // auto corrected_read = call_corrector.operator()<bio::FastaRecord<false>>();
  } else {
    corrected_read = call_corrector.operator()<bio::FastqRecord<false>>();
  }  
  auto end = std::chrono::steady_clock::now();

  spdlog::info(
      "Corrected {} reads in {} ms", corrected_read.size(),
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count());

  // auto fout = std::ofstream(output_path, std::ios::out);
  // for (auto& r : corrected_read) {
  //   fout << r << std::endl;
  // }
}