#include <filesystem>
#include <fstream>

#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>
#include <spoa/spoa.hpp>
#include <omp.h>

#include "algo/correction.hpp"
#include "utility/all.hpp"

namespace bpo = boost::program_options;
namespace fs = std::filesystem;

auto check_argument(const bpo::variables_map& vmap) {
  auto reads_path = vmap["raw_read"].as<fs::path>();
  auto overlap_path = vmap["overlap_info"].as<fs::path>();
  auto output_path = vmap["output_path"].as<fs::path>();

  /* check file format */
  check_file_format("raw_read", reads_path,
                    FILE_FORMAT::FASTA | FILE_FORMAT::FASTQ);
  check_file_format("overlap_info", overlap_path, FILE_FORMAT::PAF);
  auto fout = std::ofstream(output_path, std::ios::out);
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
  spdlog::set_level(spdlog::level::trace);
  /**
   * 1. raw reads
   * 2. overlap info
   * 3. output path
   * 4. sequencing platform
   */

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
        "platform,p", bpo::value<std::string>()->default_value("PacBio"),
        "the sequencing platform of raw reads")(
          "thread,t", bpo::value<int>()->default_value(1),
          "the number of threads to use"
        );
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

  auto raw_reads = std::vector<bio::FastaRecord<false>>{};
  if (parse_file_format(reads_path) == FILE_FORMAT::FASTA) {
    raw_reads = read_records<bio::FastaRecord<false>>(reads_path);
  } else {
    for (auto& r : read_records<bio::FastqRecord<false>>(reads_path)) {
      auto record = bio::FastaRecord<false>{.name = std::move(r.name),
                                            .seq = std::move(r.seq)};
      raw_reads.emplace_back(std::move(record));
    }
  }
  auto overlap = read_records<bio::PafRecord>(overlap_path);
  auto thread = vmap["thread"].as<int>();
  omp_set_num_threads(thread);
  

  auto start = std::chrono::steady_clock::now();
  auto correcter = FragmentedReadCorrector(std::move(raw_reads), std::move(overlap), platform, thread);
  auto corrected_read = correcter.correct();
  auto end = std::chrono::steady_clock::now();

  spdlog::info(
      "Corrected {} reads in {} ms", corrected_read.size(),
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count());

  auto fout = std::ofstream(output_path, std::ios::out);
  for (auto& r : corrected_read) {
    fout << r;
  }
}