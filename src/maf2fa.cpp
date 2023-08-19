#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "utility/all.hpp"

namespace bio = biovoltron;
namespace fs = std::filesystem;
namespace bpo = boost::program_options;

auto check_arguments(bpo::variables_map& vmap) {
  bpo::notify(vmap);
  auto fa = vmap["fa"].as<fs::path>();
  auto fout = std::ofstream(fa, std::ios::out);
  if (!fout.is_open()) {
    spdlog::error("Cannot open output file: {}", fa.string());
    exit(1);
  }
  auto thread = vmap["thread"].as<int>();
  if (thread < 1) {
    spdlog::error("Thread must be greater than 0");
    exit(1);
  }
}

int main(int argc, char* argv[]) {

  bpo::options_description opts{};
  bpo::variables_map vmap;
  try {
    opts.add_options()("help,h", "Show help message")
      ("ref,r",
        bpo::value<fs::path>()->required()->notifier(make_path_checker("ref")),
        "reference genome")
      ("maf,m",
        bpo::value<fs::path>()->required()->notifier(make_path_checker("maf")),
        "path to maf file")
      ("fa,o", bpo::value<fs::path>()->required(),
        "path to output fasta file")
      ("thread,t", bpo::value<int>()->default_value(1), "number of threads");
    bpo::store(bpo::parse_command_line(argc, argv, opts), vmap);
    check_arguments(vmap);
  } catch (const std::exception& e) {
    spdlog::error(e.what());
    exit(1);
  }

  auto ref_path = vmap["ref"].as<fs::path>();
  auto maf_path = vmap["maf"].as<fs::path>();
  auto fa_path = vmap["fa"].as<fs::path>();

  auto ref = read_records<bio::FastaRecord<false>>(ref_path);
  auto [maf_header, maf_records] = read_maf(maf_path);
  auto perfect_reads = std::vector<bio::FastaRecord<false>>{};
  for (auto& maf : maf_records) {
    auto name = maf.aln[0].name;
    auto st = maf.ref.start;
    auto len = maf.ref.len;

    auto perfect_seq = ref[0].seq.substr(st, len);
    auto perfect_read = bio::FastaRecord<false>{
      .name = name,
      .seq = ref[0].seq.substr(st, len)
    };
    perfect_reads.emplace_back(std::move(perfect_read));
  }

  auto fout = std::ofstream(fa_path, std::ios::out);
  for (const auto& read : perfect_reads) {
    fout << read << std::endl;
  }
}