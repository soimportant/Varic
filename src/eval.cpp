
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

#include <biovoltron/file_io/all.hpp>
#include <boost/program_options.hpp>
#include <omp.h>
#include <spoa/spoa.hpp>
#include <spdlog/spdlog.h>
#include "edlib.h"
#include "utility/all.hpp"

namespace bpo = boost::program_options;
namespace fs = std::filesystem;
namespace bio = biovoltron;

/**
 * @brief check arguments that passed to this program
 * @param vmap boost::program_options::variables_map
 * @return None
 */
auto check_arguments(bpo::variables_map& vmap) {
  bpo::notify(vmap);
  auto ploidy = vmap["ploidy"].as<std::size_t>();

  // the amount of path must match the ploidy and the format should be correct
  auto check_amount_and_format = [&](const std::string& opt_name,
                                     const std::vector<fs::path>& paths,
                                     const int accept_format) {
    if (paths.size() != ploidy) {
      spdlog::error("Invalid number of arguments, expected {} arguments",
                    ploidy);
      throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                  opt_name);
    }
    for (const auto& path : paths) {
      if (!fs::exists(path)) {
        spdlog::error("File {} does not exist", path.string());
        throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                    opt_name);
      }
      check_file_format(opt_name, path, accept_format);
    }
  };

  check_file_format("raw_read", vmap["raw_read"].as<fs::path>(),
                    FILE_FORMAT::FASTQ);
  for (auto& read : vmap["corrected_read"].as<std::vector<fs::path>>()) {
    check_file_format("corrected_read", read, FILE_FORMAT::FASTA);
  }
  check_amount_and_format("reference",
                          vmap["reference"].as<std::vector<fs::path>>(),
                          FILE_FORMAT::FASTA);
  check_amount_and_format("maf", vmap["maf"].as<std::vector<fs::path>>(),
                          FILE_FORMAT::MAF);
  check_amount_and_format(
      "snp_vcf", vmap["snp_vcf"].as<std::vector<fs::path>>(), FILE_FORMAT::VCF);
  check_amount_and_format("indel_vcf",
                          vmap["indel_vcf"].as<std::vector<fs::path>>(),
                          FILE_FORMAT::VCF);
  auto platform = vmap["platform"].as<std::string>();
  if (platform != "ONT" && platform != "PacBio") {
    spdlog::error("Invalid sequencing platform, only accept ONT or PacBio");
    throw bpo::validation_error(bpo::validation_error::invalid_option_value,
                                "platform");
  }
  auto output_dir = vmap["output_dir"].as<fs::path>();
  if (!fs::exists(output_dir)) {
    spdlog::info("Output directory {} does not exist", output_dir.string());
    spdlog::info("Create output directory {}", output_dir.string());
    fs::create_directories(output_dir);
  }
}

auto eval_cigar(const bio::Cigar& cigar) {
  constexpr auto eval_len = 2000ul;
  auto clen = 0ul, ac = 0ul;
  for (const auto& [len, op] : cigar | std::views::reverse) {
    if (op == '=') {
      ac += std::min(eval_len - clen, (std::size_t) len);
    }
    clen += len;
    if (clen >= eval_len) {
      break;
    }
  }
  return ac / (double) std::min(eval_len, clen);
}

/**
 * @brief evaluate the corrected read on one haplotype
 *
 * @param haplotype haplotype contains all data
 * @param output_path the output path of evaluation result
 * @return auto
 */
auto eval(const Haplotype& haplotype, const fs::path& output_path) {
  auto ref = std::string_view(haplotype.reference.seq);

  /**
   * For each read, scan the variant that covered by the read. Check whether
   * the variant is fixed correctly or not.
   *
   * Also, evaluate the correction result in front side and rear side of read
   */

  std::vector<CorrectedResult> results(haplotype.corrected_reads.size());
  spdlog::info("Corrected read size = {}", results.size());

  auto eval_one_read = [&](std::size_t i) {
    auto corrected_read = haplotype.corrected_reads[i];
    auto it = std::ranges::upper_bound(haplotype.raw_reads, corrected_read.name, {},
                                       &bio::FastqRecord<false>::name);
    assert(it != haplotype.raw_reads.begin());
    it = std::ranges::prev(it);
    const auto& raw_read = *it;
    if (corrected_read.name.find(raw_read.name) == std::string::npos) {
      spdlog::debug("raw_read = {}, corrected_read = {}", raw_read.name,
                    corrected_read.name);
    }
    auto maf = get_maf_record(haplotype.maf.second, raw_read.name);
    auto ref_start = maf.ref.start;
    auto ref_end = maf.ref.start + maf.ref.len - 1;
    auto perfect_read = ref.substr(maf.ref.start, maf.ref.len);

    auto res = CorrectedResult(maf.aln[0].strand == '+'
                                   ? std::string(perfect_read)
                                   : bio::Codec::rev_comp(perfect_read),
                               raw_read, corrected_read);
    const auto& snp_vcf_record = haplotype.snp_vcf.second;
    const auto& indel_vcf_record = haplotype.indel_vcf.second;
    res.correction_eval();
    res.variant_eval(get_vcf_range(snp_vcf_record, ref_start, ref_end), maf);
    res.variant_eval(get_vcf_range(indel_vcf_record, ref_start, ref_end), maf);
    results[i] = std::move(res);
  };

  auto tp = bio::make_threadpool(omp_get_max_threads());
  std::vector<std::future<void>> futures;
  for (auto i : std::views::iota(0u, results.size())) {
    auto [_, res] = tp.submit(eval_one_read, i);
    futures.emplace_back(std::move(res));
  }
  for (auto& f : futures) {
    f.get();
  }
  save_result(results, output_path);
}

// #pragma omp parallel for reduction(+ : v_mismatch, v_out_of_range, v_total,
//                                        total_acc, total_base)
//   for (const auto& read : haplotype.corrected_reads) {
//     auto m = get_maf_record(haplotype.maf.second, read.name);
//     auto ref_seq = ref.substr(m.ref.start, m.ref.len);
//     auto read_seq =
//         m.aln[0].strand == '+' ? read.seq : bio::Codec::rev_comp(read.seq);

//     auto cigar = align(read_seq, ref_seq);
//     for (const auto [len, op] : cigar) {
//       if (op == 'X') {
//         e_mismatch[len] += 1;
//       } else if (op == 'I') {
//         e_insertion[len] += 1;
//       } else if (op == 'D') {
//         e_deletion[len] += 1;
//       }
//     }
//     // spdlog::debug("cigar = {}", std::string(cigar));

//     total_acc += eval_cigar(cigar);
//     total_base += m.ref.len;

//     /* calculate the offset based on alignment result(cigar string) */
//     auto move_offset = cal_offset(cigar);

//     for (auto v : get_vcf_range(haplotype.snp_vcf.second, m.ref.start,
//                                 m.ref.start + m.ref.len - 1)) {
//       /* get variant start and end related on haplotype coordinate system
//       */
//       /* change all coordinate to 0-based */
//       auto vinfo = parse_vcf_info(v);
//       auto st = std::stoul(std::string(vinfo["sim_start"])) - 1 -
//       m.ref.start; auto ed = std::stoul(std::string(vinfo["sim_end"])) - 1
//       - m.ref.start; st = transform_coordinate(move_offset, st); ed =
//       transform_coordinate(move_offset, ed); v_total += 1; if (st >=
//       read_seq.size()) {
//         v_out_of_range += 1;
//         break;
//       }
//       auto alt_in_cread = read_seq.substr(st, ed - st + 1);
//       if (v.alt != alt_in_cread) {
//         v_mismatch += 1;
//       }
//     }
//   }

//   // spdlog::info("Total {} bases", total_base);
//   auto sum = 0ull;
//   auto print_error = [&](const std::map<std::size_t, std::atomic_int>& e) {
//     const int k = 10;
//     std::vector<int> v(k + 1);
//     for (const auto& [len, cnt] : e) {
//       sum += 1ull * cnt * len;
//       if (len >= k) {
//         v[k] += cnt;
//       } else {
//         v[len] += cnt;
//       }
//     }
//     for (int i = 1; i <= k; i++) {
//       spdlog::info("{}\t{}", i, v[i]);
//     }
//   };
//   spdlog::info("Mismatch");
//   print_error(e_mismatch);
//   spdlog::info("Insertion");
//   print_error(e_insertion);
//   spdlog::info("Deletion");
//   print_error(e_deletion);
//   spdlog::info("Error rate = {}", (double) sum / total_base);

//   spdlog::info("Total {} corrected reads that belong to this haplotype",
//                haplotype.corrected_reads.size());
//   spdlog::info("Covered {} variants", v_total);
//   spdlog::info(
//       "Mismatch {} variants(corrected result not match alt in vcf file)",
//       v_mismatch);
//   spdlog::info("Out of range {} variants(at tail that unavilable to
//   align)",
//                v_out_of_range);
//   spdlog::info("Accuracy at len 10 = {}",
//                (double) total_acc / haplotype.corrected_reads.size());

int main(int argc, char* argv[]) {
  spdlog::set_level(spdlog::level::debug);
  /**
   * input:
   * 1. corrected read
   * 2. original read
   * 3. reference genome
   * 4. vcf of reference genome
   *
   *
   * TODO:
   * 1. calculate the number of different variant type
   * 2. for each variant type, calculate its precision, recall, f1-score
   *
   */

  // TODO: add threads option
  bpo::options_description opts{};
  bpo::variables_map vmap;
  try {
    opts.add_options()("help,h", "Show help message")(
        "ploidy,n", bpo::value<std::size_t>()->required())(
        "raw_read,a", bpo::value<fs::path>()->required()->notifier(
                          make_path_checker("raw_read")))(
        "corrected_read,c",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "reference,r",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "maf", bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "snp_vcf",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "indel_vcf",
        bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "platform,p", bpo::value<std::string>()->required())(
        "threads,t", bpo::value<std::size_t>()->default_value(1))(
        "output_dir,o", bpo::value<fs::path>()->required());
    bpo::store(bpo::parse_command_line(argc, argv, opts), vmap);
    check_arguments(vmap);
  } catch (std::exception& ex) {
    spdlog::error("{}", ex.what());
    std::cerr << opts << std::endl;
    exit(-1);
  }

  auto ploidy = vmap["ploidy"].as<std::size_t>();
  auto haplotypes = std::vector<Haplotype>(ploidy);
  auto raw_reads = std::vector<bio::FastqRecord<false>>{};
  auto threads = vmap["threads"].as<std::size_t>();
  auto output_dir = vmap["output_dir"].as<fs::path>();

  omp_set_num_threads(threads);

  /* --- read data --- */
  for (auto i : std::views::iota(0u, ploidy)) {
    auto ref_path = vmap["reference"].as<std::vector<fs::path>>()[i];
    auto maf_path = vmap["maf"].as<std::vector<fs::path>>()[i];
    auto snp_path = vmap["snp_vcf"].as<std::vector<fs::path>>()[i];
    auto indel_path = vmap["indel_vcf"].as<std::vector<fs::path>>()[i];

    haplotypes[i].reference =
        read_records<bio::FastaRecord<false>>(ref_path)[0];
    haplotypes[i].maf = read_maf(maf_path);
    haplotypes[i].snp_vcf = read_vcf(snp_path);
    haplotypes[i].indel_vcf = read_vcf(indel_path);
  }
  /* --- data preprocessing --- */
  /* for vechat, it will add "rr" suffix to corrected read name, holy moly!!! */
  // std::ranges::for_each(corrected_reads, [](auto &r) {
  //   r.name = r.name.substr(0, r.name.size() - 2);
  // });
  // std::ranges::sort(raw_reads, {}, &bio::FastqRecord<false>::name);
  // std::ranges::sort(corrected_reads, {}, &bio::FastaRecord<false>::name);
  std::ranges::for_each(haplotypes, [](auto& h) {
    std::ranges::sort(h.maf.second, [](const auto& lhs, const auto& rhs) {
      return lhs.aln[0].name < rhs.aln[0].name;
    });
  });

  auto raw_read_path = vmap["raw_read"].as<fs::path>();
  raw_reads = read_records<bio::FastqRecord<false>>(raw_read_path);
  std::sort(
      std::execution::par, raw_reads.begin(), raw_reads.end(),
      [](const auto& lhs, const auto& rhs) { return lhs.name < rhs.name; });
  for (auto p : std::views::iota(0u, ploidy)) {
    auto read_prefix = std::to_string(p);
    std::ranges::for_each(raw_reads, [&](auto& r) {
      if (r.name.starts_with(read_prefix)) {
        haplotypes[p].raw_reads.emplace_back(std::move(r));
      }
    });
  }
  for (const auto& path : vmap["corrected_read"].as<std::vector<fs::path>>()) {
    spdlog::info("Evalation on {}", path.string());
    auto corrected_reads = read_records<bio::FastaRecord<false>>(path);
    /* seperate raw and corrected reads according to prefix of read name */
    for (auto p : std::views::iota(0u, ploidy)) {
      auto read_prefix = fmt::format("{}", p);
      haplotypes[p].corrected_reads.clear();
      std::ranges::for_each(corrected_reads, [&](auto& r) {
        if (r.name.starts_with(read_prefix)) {
          r.name = r.name.substr(0, r.name.size() - 2); // ! only for vechat
          haplotypes[p].corrected_reads.emplace_back(std::move(r));
        }
      });
    }
    std::vector<CorrectedResult> results(ploidy);
    for (const auto p : std::views::iota(0u, ploidy)) {
      auto output_path = output_dir / fmt::format("h{}_result.csv", p);
      spdlog::info("Evalation on haplotype {}, write result into {}", p, output_path.string());
      eval(haplotypes[p], output_path);
    }
  }
}