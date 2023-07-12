#include <spdlog/spdlog.h>

#include <biovoltron/file_io/all.hpp>
#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>
#include <spoa/spoa.hpp>

#include "edlib.h"
#include "utility/all.hpp"

namespace bpo = boost::program_options;
namespace fs = std::filesystem;
namespace bio = biovoltron;

auto check_arguments(bpo::variables_map &vmap) {
  bpo::notify(vmap);
  auto ploidy = vmap["ploidy"].as<std::size_t>();
  // the amount of path must match the ploidy and the format must acceptable
  auto check_amount_and_format = [&](const std::string &opt_name,
                                     const std::vector<fs::path> &paths,
                                     const int accept_format) {
    if (paths.size() != ploidy) {
      throw bpo::validation_error(
          bpo::validation_error::invalid_option_value, opt_name,
          fmt::format("Invalid number of arguments, expected {} arguments",
                      ploidy));
    }
    for (const auto &path : paths) {
      if (!fs::exists(path)) {
        throw bpo::validation_error(
            bpo::validation_error::invalid_option_value, opt_name,
            fmt::format("File {} does not exist", path.string()));
      }
      check_file_format(opt_name, path, accept_format);
    }
  };

  check_file_format("raw_read", vmap["raw_read"].as<fs::path>(),
                    FILE_FORMAT::FASTQ);
  check_file_format("corrected_read", vmap["corrected_read"].as<fs::path>(),
                    FILE_FORMAT::FASTA);
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
}

auto get_cigar(std::string_view ref, std::string_view read) {
  auto result = edlibAlign(
      ref.data(), ref.size(), read.data(), read.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  auto cigar = std::unique_ptr<char>(edlibAlignmentToCigar(
      result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED));
  edlibFreeAlignResult(result);
  return bio::Cigar(cigar.get());
};

/**
 * @brief Use for calculate the relative offset between reference and read, the
 * cigar string must be the alignment result from **read to reference**.
 *
 * @param cigar Cigar string between reference and read
 * @return std::vector<std::pair<std::size_t, int>>
 */
auto cal_offset(const bio::Cigar &cigar, bool rev = false) {
  std::vector<std::pair<std::size_t, std::int64_t>> offset;
  offset.emplace_back(0u, 0);
  for (auto i = 0u, haplotype_pos = 0u; i < cigar.size(); i += 1) {
    auto dx = 0;
    if (cigar[i].op == 'I') {
      dx = cigar[i].size;
    } else {
      if (cigar[i].op == 'D') {
        dx = -cigar[i].size;
      }
      haplotype_pos += cigar[i].size;
    }
    if (dx) {
      dx += offset.back().second;
      offset.emplace_back(haplotype_pos, dx);
    }
  }
  if (rev) {
    for (auto& [p, k] : offset) {
      p += k;
      k = -k;
    }
  }
  return offset;
}

/**
 * @brief Calculate the coordinate transformation from reference genome to
 * haplotype by using vcf contains only indels.
 * @param indels vcf files contains only indels
 * @param rev
 * @return
 */
auto cal_offset(const std::vector<bio::VcfRecord> &indels, bool rev = false) {
  std::vector<std::pair<std::size_t, std::int64_t>> offset;
  offset.emplace_back(0u, 0);
  /**
   * NC_000913.3	8216	.	C	CGT	.	. variant_type=INDEL;
   * NC_000913.3	9190	.	A	ACA	.	. variant_type=INDEL;
   */
  // for (auto ref_pos = 0u; auto &v : indels) {
  //   auto dx = std::ranges::ssize(v.alt) - std::ranges::ssize(v.ref);
  //   if (rev) {
  //     auto hpos = v.pos - offset.back().second;
  //     offset.emplace_back(hpos, offset.back().second + dx);
  //   } else {
  //     offset.emplace_back(v.pos + dx, offset.back().second - dx);
  //   }
  // }
  for (auto &v : indels) {
    // ssize -> signed size()
    auto dx = std::ranges::ssize(v.alt) - std::ranges::ssize(v.ref);
    offset.emplace_back(v.pos + dx, offset.back().second - dx);
  }
  if (rev) {
    for (auto& [p, k] : offset) {
      p += k;
      k = -k;
    }
  }
  return offset;
}

auto transform_coordinate(
    const std::vector<std::pair<std::size_t, std::int64_t>> &offset,
    std::size_t old_pos) {
  auto it = std::ranges::upper_bound(
      offset, old_pos, {}, &std::pair<std::size_t, std::int64_t>::first);
  assert(it != offset.begin());
  return old_pos + std::prev(it)->second;
}

auto transform_pair_coordinate(
    const std::vector<std::pair<std::size_t, std::int64_t>> &offset,
    std::size_t st, std::size_t ed) {
  auto it_st = std::ranges::upper_bound(
      offset, st, {}, &std::pair<std::size_t, std::int64_t>::first);
  auto it_ed = std::ranges::lower_bound(
      offset, ed, {}, &std::pair<std::size_t, std::int64_t>::first);
  assert(it_st != offset.begin() && it_ed != offset.begin());
  return std::make_pair(st + std::prev(it_st)->second, 
                        ed + std::prev(it_ed)->second);
}

auto eval_cigar(const bio::Cigar &cigar) {
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

class Haplotype {
public:
  Haplotype() = default;
  bio::FastaRecord<false> reference;
  std::vector<bio::FastqRecord<false>> raw_reads;
  std::vector<bio::FastaRecord<false>> corrected_reads;
  std::pair<bio::VcfHeader, std::vector<bio::VcfRecord>> snp_vcf;
  std::pair<bio::VcfHeader, std::vector<bio::VcfRecord>> indel_vcf;
  std::pair<bio::MafHeader, std::vector<bio::MafRecord>> maf;
};

auto eval(const Haplotype &haplotype) {
  auto ref = std::string_view(haplotype.reference.seq);

  /* evaluate variant */
  auto v_mismatch = 0u;
  auto v_out_of_range = 0u;
  auto v_total = 0u;

  /* evaluate correction result */
  auto total_acc = (double) 0;
  auto e_mismatch = std::map<std::size_t, std::atomic_int>{};
  auto e_insertion = std::map<std::size_t, std::atomic_int>{};
  auto e_deletion = std::map<std::size_t, std::atomic_int>{};

  #pragma omp parallel for reduction(+:v_mismatch, v_out_of_range, v_total, total_acc)
  for (const auto &read : haplotype.corrected_reads) {
    auto m = get_maf_record(haplotype.maf.second, read.name);
    auto ref_seq = ref.substr(m.ref.start, m.ref.len);
    auto read_seq =
        m.aln[0].strand == '+' ? read.seq : bio::Codec::rev_comp(read.seq);

    auto cigar = get_cigar(read_seq, ref_seq);
    for (const auto [len, op] : cigar) {
      if (op == 'X') {
        e_mismatch[len] += 1;
      } else if (op == 'I') {
        e_insertion[len] += 1;
      } else if (op == 'D') {
        e_deletion[len] += 1;
      }
    }
    total_acc += eval_cigar(cigar);

    /* calculate the offset based on alignment result(cigar string) */
    auto move_offset = cal_offset(cigar);

    for (auto v : get_vcf_range(haplotype.snp_vcf.second, m.ref.start,
                                m.ref.start + m.ref.len - 1)) {
      /* get variant start and end related on haplotype coordinate system */
      /* change all coordinate to 0-based */
      auto vinfo = parse_vcf_info(v);
      auto st = std::stoul(std::string(vinfo["sim_start"])) - 1 - m.ref.start;
      auto ed = std::stoul(std::string(vinfo["sim_end"])) - 1 - m.ref.start;
      st = transform_coordinate(move_offset, st);
      ed = transform_coordinate(move_offset, ed);
      v_total += 1;
      if (st >= read_seq.size()) {
        v_out_of_range += 1;
        break;
      }
      auto alt_in_cread = read_seq.substr(st, ed - st + 1);
      if (v.alt != alt_in_cread) {
        v_mismatch += 1;
      }
    }

    // spdlog::info("Total {} corrected reads that belong to this haplotype",
    //              haplotype.corrected_reads.size());
    // spdlog::info("Covered {} variants", v_total);
    // spdlog::info(
    //     "Mismatch {} variants(corrected result not match alt in vcf file)",
    //     v_mismatch);
    // spdlog::info("Out of range {} variants(at tail that unavilable to align)",
    //              v_out_of_range);

    // spdlog::info("Total {} bases", total_base);
    auto sum = 0ull;
    auto print_error = [&](const std::map<std::size_t, std::atomic_int> &e) {
      const int k = 10;
      std::vector<int> v(k + 1);
      for (const auto &[len, cnt] : e) {
        sum += 1ull * cnt * len;
        if (len >= k) {
          v[k] += cnt;
        } else {
          v[len] += cnt;
        }
      }
      for (int i = 1; i <= k; i++) {
        spdlog::info("{}\t{}", i, v[i]);
      }
    };
    // spdlog::info("Mismatch");
    // print_error(e_mismatch);
    // spdlog::info("Insertion");
    // print_error(e_insertion);
    // spdlog::info("Deletion");
    // print_error(e_deletion);
    // spdlog::info("Error rate = {}", (double) sum / total_base);
  }
  spdlog::info("Accuracy at len 10 = {}",
               (double) total_acc / haplotype.corrected_reads.size());
}

// auto eval(const Haplotype &h1, const Haplotype &h2 /* may have reuslt here */) {
//   auto h1_to_ref = cal_offset(h1.indel_vcf.second);
//   auto ref_to_h1 = cal_offset(h1.indel_vcf.second, true);
//   auto h2_to_ref = cal_offset(h2.indel_vcf.second);
//   auto ref_to_h2 = cal_offset(h2.indel_vcf.second, true);
//   auto h1_to_h2 = [&](std::size_t pos) {
//     pos = transform_coordinate(h1_to_ref, pos);
//     pos = transform_coordinate(ref_to_h2, pos);
//     return pos;
//   };
//   auto h2_to_h1 = [&](std::size_t pos) {
//     pos = transform_coordinate(h2_to_ref, pos);
//     pos = transform_coordinate(ref_to_h1, pos);
//     return pos;
//   };

//   auto h1_ref = std::string_view(h1.reference.seq);
//   auto h2_ref = std::string_view(h2.reference.seq);

//   auto total_variant = 0;
//   auto wrong_ins = 0;
//   auto wrong_del = 0;
//   auto wrong_snp = 0;
//   auto out_of_range = 0;

//   #pragma omp parallel for reduction(+:total_variant, wrong_ins, wrong_del, out_of_range)
//   for (auto h1_cread : h1.corrected_reads) {
//     auto m = get_maf_record(h1.maf.second, h1_cread.name);
//     auto read_seq = m.aln[0].strand == '+' ? h1_cread.seq
//                                            : bio::Codec::rev_comp(h1_cread.seq);
//     // the base of this read from h1 is a variant on h2
//     // -> I need to know the variants that covered by this read in h2
//     // -> we have the read range on haplotype1

//     // -> for each variant on h2:
//     // ->   find the corresponding base on h1
//     // ->   check whether it is same as h1
//     // ->     if match -> correct successfully
//     //          otherwise -> correct failed -> wrong haplotype

//     auto print_genome = [&](std::string_view seq, std::size_t pos) {
//       constexpr std::size_t len = 10;
//       std::size_t st = std::max(0L, (std::int64_t) pos - (std::int64_t) len);
//       std::size_t ed = std::min(seq.size(), pos + len);
//       spdlog::debug("{:10} : {:>21}", st, seq.substr(st, ed - st + 1));
//       spdlog::debug("{:10}   {:>11}", "", '^');
//     };

//     /**
//      * Method 1
//      * 1. get corresponding corrdinate of read on h1
//      * 2. transform coordinate from h1 to h2
//      * 3. for each variant on h2, reverse it's coordinate to read
//      * 4. check read subseq == h2 subseq
//      */
//     // {
//     //   auto cigar = get_cigar(h1_ref.substr(m.ref.start, m.ref.len), read_seq);
//     //   auto read_to_h1 = cal_offset(cigar);
//     //   auto h1_to_read = cal_offset(cigar, true);
//     //   auto h2_st = h1_to_h2(m.ref.start);
//     //   auto h2_ed = h1_to_h2(m.ref.start + m.ref.len - 1);

//     //   for (const auto& v : get_vcf_range(h2.snp_vcf.second, h2_st, h2_ed)) {
//     //     auto vinfo = parse_vcf_info(v);

//     //     auto v_h2_st = std::stoul(std::string(vinfo.at("sim_start"))) - 1;
//     //     auto v_h2_ed = std::stoul(std::string(vinfo.at("sim_end"))) - 1;
//     //     auto v_h1_st = h2_to_h1(v_h2_st);
//     //     auto v_h1_ed = h2_to_h1(v_h2_ed);
//     //     auto v_read_st = transform_coordinate(h1_to_read, v_h1_st - m.ref.start);
//     //     auto v_read_ed = transform_coordinate(h1_to_read, v_h1_ed - m.ref.start);

//     //     auto v_read_seq = read_seq.substr(v_read_st, v_read_ed - v_read_st + 1);
//     //     auto v_h2_seq = h2_ref.substr(v_h2_st, v_h2_ed - v_h2_st + 1);
//     //     auto v_h1_seq = h1_ref.substr(v_h1_st, v_h1_ed - v_h1_st + 1);

//     //     if (v_read_seq == v_h2_seq) {
//     //       spdlog::debug("v.pos = {}", v.pos);
//     //       spdlog::info("Correct to wrong haplotype!!!");
//     //       spdlog::info("{}", v_h2_seq);
//     //       spdlog::info("{}", v_read_seq);
//     //       spdlog::info("{}", v_h1_seq);

//     //       spdlog::debug("v_h2_st: {}", v_h2_st);
//     //       spdlog::debug("v_h1_st: {}", v_h1_st);
//     //       spdlog::debug("v_read_st: {}", v_read_st);

//     //       print_genome(h2_ref, v_h2_st);
//     //       print_genome(read_seq, v_read_st);
//     //       print_genome(h1_ref, v_h1_st);
//     //     }
//     //   }
//     // }

//     /**
//      * Method 2
//      * 1. get corresponding corrdinate of read on h1
//      * 2. transform coordinate from h1 to h2
//      * 3. get alignment result from read to h2
//      * 4. for each variant in h2, check whether it exists in read
//      */
//     {
//       auto h2_st = h1_to_h2(m.ref.start);
//       auto h2_ed = h1_to_h2(m.ref.start + m.ref.len - 1);
//       auto cigar = get_cigar(h2_ref.substr(h2_st, h2_ed - h2_st + 1), read_seq);
//       auto cigar_seg = std::vector<std::size_t>(cigar.size() + 1);
//       auto h2_to_read = cal_offset(cigar, true);

//       cigar_seg[0] = 0;
//       for (auto i = 1u; i <= cigar_seg.size(); i++) {
//         cigar_seg[i] += cigar_seg[i - 1];
//       }
//       auto get_ele_from_cigar = [&](std::size_t pos) {
//         auto it = std::ranges::lower_bound(cigar_seg, pos);
//         if (*it != pos) {
//           assert(it != cigar_seg.begin());
//           it = std::prev(it);
//         }
//         auto d = std::ranges::distance(cigar_seg.begin(), it);
//         assert(d != 0);
//         return d - 1;
//       };

//       // spdlog::debug("cigar = {}", std::string(cigar));
//       for (const auto& v : get_vcf_range(h2.indel_vcf.second, h2_st, h2_ed)) {
//         auto vinfo = parse_vcf_info(v);
//         /* remember that coordinate of vcf is 1-based */
//         auto v_h2_st = std::stoul(std::string(vinfo.at("sim_start"))) - 1;
//         auto v_h2_ed = std::stoul(std::string(vinfo.at("sim_end"))) - 1;
//         auto v_read_st = transform_coordinate(h2_to_read, v_h2_st - h2_st);
//         auto v_read_ed = transform_coordinate(h2_to_read, v_h2_ed - h2_st);

//         // auto [v_read_st, v_read_ed] = 
//         //   transform_pair_coordinate(h2_to_read, v_h2_st - h2_st, v_h2_ed - h2_st);

//         // assert(v_read_st <= v_read_ed);
//         total_variant++;
//         if (v_read_st > read_seq.size() || v_read_st > v_read_ed) {
//           out_of_range++;
//           continue;
//         }

//         auto v_read_seq = read_seq.substr(v_read_st, v_read_ed - v_read_st + 1);
//         auto v_h2_seq = h2_ref.substr(v_h2_st, v_h2_ed - v_h2_st + 1);

//         /**
//          * NC_000913.3	568985	.	TG	T	.	. variant_type=INDEL
//          * NC_000913.3	568986	.	G	GACA	.	. variant_type=INDEL
//          */

//         // {
//         //   auto ref_path = "/mnt/ec/ness/yolkee/thesis/data/ref/Ecoli/K12/"
//         //                   "GCF_000005845.2/ncbi_dataset/data/GCF_000005845.2/"
//         //                   "GCF_000005845.2_ASM584v2_genomic.fna";
//         //   auto ref_fa = read_records<bio::FastaRecord<false>>(ref_path);
//         //   auto ref = std::string_view(ref_fa[0].seq);
//         //   spdlog::debug("v_h2_st = {}", v_h2_st);
//         //   spdlog::debug("v_h2_st - h2_st = {}", v_h2_st - h2_st);
//         //   spdlog::debug("v_read_st = {}", v_read_st);
//         //   spdlog::debug("v.alt = {}, v_h2_seq = {}", v.alt, v_h2_seq);
//         //   print_genome(h2_ref, v_h2_st);
//         //   auto ref_st = std::stoul(std::string(vinfo.at("ref_start"))) - 1;
//         //   auto ref_end = std::stoul(std::string(vinfo.at("ref_end"))) - 1;
//         //   print_genome(ref, ref_st);
//         //   print_genome(read_seq, v_read_st);
//         // }

//         /**
//          * If the variant type is deletion, and the alignment result of alt in
//          * this variant and following base should be match, then it can be 
//          * consider as correct to wrong haplotype.
//          * 
//          * e.g. 
//          * v.ref = TAT
//          * v.alt = T
//          * 
//          * ref   = ATGTATTGA
//          * h     = ATGT--TGA
//          * read  = ATGTA-TGA
//          * cigar = 4M1I3M
//          * 
//          * the variant is on 4M
//          */
        
//         // if (v_read_seq.size() > 20) {
//         //   spdlog::debug("cigar = {}", std::string(cigar));
//         //   spdlog::debug("v_read_seq = {}", v_read_seq);
//         //   spdlog::debug("v_h2_seq = {}", v_h2_seq);
//         //   spdlog::debug("v_read = ({}, {})", v_read_st, v_read_ed);
//         //   spdlog::debug("v_h2 = ({}, {})", v_h2_st - h2_st, v_h2_ed - h2_st);
//         //   for (const auto& [pos, offset] : h2_to_read) {
//         //     spdlog::debug("h2_to_read = ({}, {})", pos, offset);
//         //   }
//         //   print_genome(h2_ref, v_h2_st);
//         //   print_genome(read_seq, v_read_st);
//         // }

//         if (v.ref.size() < v.alt.size()) {
//           /* Insertion */
//           auto sz_diff = (std::int64_t) (v.ref.size() - v.alt.size());
//           auto next_non_match = v_read_ed;
//           while (read_seq[next_non_match] == v.alt.back()) {
//             next_non_match += 1;
//           }
//           auto c = get_ele_from_cigar(next_non_match);
//           assert(c != cigar.size());
//           if (cigar[c].op != 'D' || cigar[c].size != sz_diff) {
//             // spdlog::info("Correct to wrong haplotype!!!");
//             wrong_ins += 1;
//           }
//         } else if (v.ref.size() > v.alt.size()) {
//           /* Deletion */
//           if (read_seq.substr(v_read_st, v.ref.size()) != v.ref) {
//             spdlog::debug("v.ref <-> read_seq = ({}, {})", v.ref, read_seq.substr(v_read_st, v.ref.size()));
//             spdlog::debug("v.alt <-> read_seq = ({}, {})", v.alt, read_seq.substr(v_read_st, v.alt.size()));
//             // spdlog::info("Correct to wrong haplotype!!!");
//             wrong_del += 1;
//           }
//         } else {
//           if (read_seq == v.alt) {
//             wrong_snp += 1;
//           }
//         }

//         // if (v.ref.size() < v.alt.size()) {
//         //   /* SNP or insertion */

//         //   static int cnt = 0;
//         //   if (cnt < 20) {
//         //     spdlog::debug("[INS] v.ref = {}, v.alt = {}", v.ref, v.alt);
//         //     cnt++;
//         //   } 

//         //   // if (v_read_seq == v_h2_seq) {
//         //   //   // spdlog::info("Correct to wrong haplotype!!!");
//         //   //   wrong_ins += 1;
//         //   // }
//         // } else {
//         //   /* Deletion */
//         //   /* if match v.ref, it can consider as correct successfully */
//         //   /* if not match */
          
//         //   /* Insertion */
//         //   /* */
//         // }

//         // if (v_read_seq == v_h2_seq) {
//         //   // spdlog::info("Correct to wrong haplotype!!!");
//         //   if (v.ref.size() < v.alt.size()) {
//         //     wrong_ins += 1;
//         //     if (wrong_ins < 20) {
//         //       // spdlog::debug("v.ref = {}, v.alt = {}", v.ref, v.alt);
//         //       // print_genome(h2_ref, v_h2_st);
//         //       // print_genome(read_seq, v_read_st);
//         //     }
//         //   } else {
//         //     // wrong_del += 1;
//         //     // if (wrong_del < 20) {
//         //     //   spdlog::debug("v.ref = {}, v.alt = {}", v.ref, v.alt);
//         //     //   print_genome(h2_ref, v_h2_st);
//         //     //   print_genome(read_seq, v_read_st);
//         //     // }
//         //   }
//         // }
//       }
//     }
//   }
//   spdlog::info("total_variant = {}", total_variant);
//   spdlog::info("out_of_range = {}", out_of_range);
//   spdlog::info("wrong_ins = {}", wrong_ins);
//   spdlog::info("wrong_del = {}", wrong_del);

// }

auto eval(const Haplotype &h1, const Haplotype& h2) {
  /**
   * 1. create chain file from h1 to h2
   * 2. align h1.corrected read to h2
   * 3. transform alignment result from h1 coordinate to h2
   */

  
  
}

int main(int argc, char *argv[]) {
  spdlog::set_level(spdlog::level::trace);
  omp_set_num_threads(7);
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


  bpo::options_description opts{};
  bpo::variables_map vmap;
  try {
    opts.add_options()("help,h", "Show help message")(
        "ploidy,n", bpo::value<std::size_t>()->required())(
        "raw_read,o", bpo::value<fs::path>()->required()->notifier(
                          make_path_checker("raw_read")))(
        "corrected_read,c", bpo::value<fs::path>()->required()->notifier(
                                make_path_checker("corrected_read")))(
        "reference,r", bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "maf", bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "snp_vcf", bpo::value<std::vector<fs::path>>()->required()->multitoken())(
        "indel_vcf", bpo::value<std::vector<fs::path>>()->required()->multitoken());
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
  auto corrected_reads = std::vector<bio::FastaRecord<false>>{};

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
  raw_reads =
      read_records<bio::FastqRecord<false>>(vmap["raw_read"].as<fs::path>());
  corrected_reads =
      read_records<bio::FastaRecord<false>>(vmap["corrected_read"].as<fs::path>());

  /* --- data preprocessing --- */
  /* for vechat, it will add "rr" suffix to corrected read name, holy moly!!! */
  std::ranges::for_each(corrected_reads, [](auto &r) {
    r.name = r.name.substr(0, r.name.size() - 2);
  });
  // std::ranges::sort(raw_reads, {}, &bio::FastqRecord<false>::name);
  // std::ranges::sort(corrected_reads, {}, &bio::FastaRecord<false>::name);
  std::ranges::for_each(haplotypes, [](auto &h) {
    std::ranges::sort(h.maf.second, [](const auto &lhs, const auto &rhs) {
      return lhs.aln[0].name < rhs.aln[0].name;
    });
  });

  /* seperate raw and corrected reads according to prefix of read name */
  for (auto p : std::views::iota(0u, ploidy)) {
    auto read_prefix = fmt::format("{}", p);
    std::ranges::for_each(raw_reads, [&](auto &r) {
      if (r.name.starts_with(read_prefix)) {
        haplotypes[p].raw_reads.emplace_back(std::move(r));
      }
    });
    std::ranges::for_each(corrected_reads, [&](auto &r) {
      if (r.name.starts_with(read_prefix)) {
        haplotypes[p].corrected_reads.emplace_back(std::move(r));
      }
    });
  }

  // for (const auto p : std::views::iota(0u, ploidy)) {
  //   eval(haplotypes[p]);
  // }
  
  eval(haplotypes[1], haplotypes[0]);

  // DataSet dataset(reference, raw_reads, corrected_reads);
  // spdlog::info("Start evaluating SNP");
  // dataset.eval("0", snp[0].first, snp[0].second, maf[0].second);
  // spdlog::info("Start evaluating INDEL");
  // dataset.eval("0", indel[0].first, indel[0].second, maf[0].second);
}