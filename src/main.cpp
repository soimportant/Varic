#include <chrono>
#include <filesystem>
#include <fstream>

#include <omp.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <boost/program_options.hpp>
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
      spdlog::error("{} does not exist", p.string());
      throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                  opt_name);
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
  if (platform.find("PacBio") == std::string::npos &&
      platform.find("ONT") == std::string::npos) {
    spdlog::error("Platform must be PacBio or ONT");
    exit(1);
  }
  auto thread = vmap["thread"].as<int>();
  if (thread < 1) {
    spdlog::error("Thread must be greater than 0");
    exit(1);
  }
  auto match = vmap["match"].as<int>();
  auto mismatch = vmap["mismatch"].as<int>();
  auto gap = vmap["gap"].as<int>();
  auto extend = vmap["extend"].as<int>();
  if (match < 1) {
    spdlog::error("Match must be greater than 0");
    exit(1);
  }
  if (mismatch > 0) {
    spdlog::error("Mismatch must be less than 0");
    exit(1);
  }
  if (gap > 0) {
    spdlog::error("Gap must be less than 0");
    exit(1);
  }
  if (extend > 0) {
    spdlog::error("Extend must be less than 0");
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
    opts.add_options()("help,h", "Show help message")
        // raw reads
        ("raw_read,r",
         bpo::value<fs::path>()->required()->notifier(
             make_path_checker("raw_read")),
         "The path to TGS raw reads")
        // overlap information
        ("overlap_info,c",
         bpo::value<fs::path>()->required()->notifier(
             make_path_checker("overlap_info")),
         "The path to the overlap information of raw reads(.paf)")
        // output path of corrected read
        ("output_path,o", bpo::value<fs::path>()->required(),
         "The output path of corrected read")
        // sequencing platform
        ("platform,p", bpo::value<std::string>(),
         "The sequencing platform of raw reads, must be \"ONT\" or \"PacBio\"")
        // Max used threads
        ("thread,t", bpo::value<int>()->default_value(1),
         "The maximum number of threads")
        // match score
        ("match", bpo::value<int>()->default_value(5),
         "The match score in alignment")
        // mismatch score
        ("mismatch", bpo::value<int>()->default_value(-4),
         "The mismatch score in alignment")
        // gap score
        ("gap", bpo::value<int>()->default_value(-8),
         "The gap score in alignment")
        // gap extension score
        ("extend", bpo::value<int>()->default_value(-6),
         "The gap extension score in alignment")
        // graph pruning threshold
        ("prune", bpo::value<double>()->default_value(0.95),
         "The prune threshold for the variation graph, this value should set "
         "between 0 and 1. The higher the value, the pruned graph will be more "
         "aggressive. If your data is noisy, which means has more haplotypes "
         "inside it, you should set this value higher")
        // maximum depth inside a window
        ("depth,d", bpo::value<int>()->default_value(-1),
         "The maximum sequence to be used when building variation graph of a "
         "window, -1 means take all sequences. This will be used to reduce the "
         "run-time of the program, but will effect the accuracy of the result.")
        // seed for random number generator
        ("seed", bpo::value<std::int64_t>(),
         "The seed for random number generator")(
            "quiet,q", bpo::bool_switch()->default_value(false),
            "disable all log")
        // debug flag
        ("debug", bpo::bool_switch()->default_value(false),
         "enable debug mode (print verbose log to stderr)");
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
  auto thread = vmap["thread"].as<int>();

  auto match = vmap["match"].as<int>();
  auto mismatch = vmap["mismatch"].as<int>();
  auto gap = vmap["gap"].as<int>();
  auto extend = vmap["extend"].as<int>();

  auto depth = vmap["depth"].as<int>();
  auto prune_ratio = vmap["prune"].as<double>();

  auto seed = vmap.contains("seed")
                  ? vmap["seed"].as<std::int64_t>()
                  : std::chrono::system_clock::now()
                        .time_since_epoch()
                        .count();
  auto quiet = vmap["quiet"].as<bool>();

  if (debug_mode) {
    auto logger = spdlog::stderr_color_mt("logger");
    spdlog::set_default_logger(logger);
    spdlog::set_level(spdlog::level::debug);
  } else {
    spdlog::set_level(spdlog::level::info);
  }
  omp_set_num_threads(thread);
  spdlog::info("Write output to {}", output_path.string());

  auto overlap = read_records<bio::PafRecord>(overlap_path);
  auto call_corrector = [&]<class R>() {
    auto raw_reads = read_records<R>(reads_path);
    auto corrector =
        FragmentedReadCorrector<R>(std::move(raw_reads), std::move(overlap),
                                   platform, thread, depth, prune_ratio, seed);
    corrector.set_alignment_params(match, mismatch, gap, extend);
    auto corrected_read = corrector.correct();
    return corrected_read;
  };

  auto corrected_read = std::vector<bio::FastaRecord<false>>{};
  if (parse_file_format(reads_path) == FILE_FORMAT::FASTA) {
    corrected_read = call_corrector.operator()<bio::FastaRecord<false>>();
  } else {
    corrected_read = call_corrector.operator()<bio::FastqRecord<false>>();
  }
  auto end = std::chrono::steady_clock::now();

  spdlog::info(
      "Corrected {} reads in {} ms", corrected_read.size(),
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count());

  auto fout = std::ofstream(output_path, std::ios::out);
  for (auto& r : corrected_read) {
    fout << r << std::endl;
  }
}