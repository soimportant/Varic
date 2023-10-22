#pragma once

#include <algorithm>
#include <array>
#include <filesystem>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "edlib.h"
#include "rapidcsv.h"
#include <biovoltron/file_io/all.hpp>
#include <spdlog/spdlog.h>

#include "thesis/utility/format/maf.hpp"
#include "thesis/utility/format/vcf.hpp"

namespace fs = std::filesystem;
namespace bio = biovoltron;

/**
 * @brief align corrected read to corresponding subsequence of reference genome
 *
 * @param ref subsequence of reference genome
 * @param read corrected read
 * @return bio::Cigar
 */
auto align(std::string_view ref, std::string_view read) {
  auto result = edlibAlign(
      read.data(), read.size(), ref.data(), ref.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  auto cigar = (edlibAlignmentToCigar(result.alignment, result.alignmentLength,
                                      EDLIB_CIGAR_EXTENDED));
  edlibFreeAlignResult(result);
  auto r = bio::Cigar(cigar);
  std::free(cigar);
  return r;
};

/**
 * @brief Use for calculate the relative offset between reference and read, the
 * cigar string must be the alignment result from **read to reference**.
 *
 * @param cigar Cigar string between reference and read
 * @return std::vector<std::pair<std::size_t, int>>
 */
auto cal_offset(const bio::Cigar& cigar, bool rev = false) {
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
 * @param rev reverse the coordinate transformation
 * @return std::vector<std::pair<std::size_t, int>>
 */
auto cal_offset(const std::vector<bio::VcfRecord>& indels, bool rev = false) {
  std::vector<std::pair<std::size_t, std::int64_t>> offset;
  offset.emplace_back(0u, 0);
  /**
   * NC_000913.3	8216	.	C	CGT	.	.
   * variant_type=INDEL; NC_000913.3	9190	.	A	ACA
   * .	. variant_type=INDEL;
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
  for (auto& v : indels) {
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

/**
 * @brief transform old coordinate to new coordinate by using offset
 *
 * @param offset offset information
 * @param old_pos old coordinate
 * @return std::size_t new coordinate
 */
auto transform_coordinate(
    const std::vector<std::pair<std::size_t, std::int64_t>>& offset,
    std::size_t old_pos) {
  auto it = std::ranges::upper_bound(
      offset, old_pos, {}, &std::pair<std::size_t, std::int64_t>::first);
  assert(it != offset.begin());
  return old_pos + std::prev(it)->second;
}

/**
 * @brief transform a pair of coordinate to new pair by using offset information
 * @param offset offset information
 * @param st start position
 * @param ed end position
 * @return auto
 */
auto transform_pair_coordinate(
    const std::vector<std::pair<std::size_t, std::int64_t>>& offset,
    std::size_t st, std::size_t ed) {
  auto it_st = std::ranges::upper_bound(
      offset, st, {}, &std::pair<std::size_t, std::int64_t>::first);
  auto it_ed = std::ranges::lower_bound(
      offset, ed, {}, &std::pair<std::size_t, std::int64_t>::first);
  assert(it_st != offset.begin() && it_ed != offset.begin());
  return std::make_pair(st + std::prev(it_st)->second,
                        ed + std::prev(it_ed)->second);
}

struct CorrectedResult {

  CorrectedResult() = default;

  CorrectedResult(const std::string& ref_seq,
                  const bio::FastqRecord<false>& raw_read,
                  const bio::FastaRecord<false>& corrected_read) {
    ref = ref_seq;
    this->raw_read = raw_read;
    this->corrected_read = corrected_read;
    cigar = align(ref, corrected_read.seq);
  }

  auto correction_eval() {
    const auto min_correct_len = 100u;
    std::size_t correct_st = 0;
    std::size_t correct_ed = cigar.size() - 1;

    while (correct_st < cigar.size()) {
      if (cigar[correct_st].op == '=' &&
          cigar[correct_st].size >= min_correct_len) {
        break;
      }
      correct_st += 1;
    }
    if (correct_st == cigar.size()) {
      has_corrected = false;
      return;
    }
    while (correct_ed > correct_st) {
      if (cigar[correct_ed].op == '=' &&
          cigar[correct_ed].size >= min_correct_len) {
        break;
      }
      correct_ed -= 1;
    }
    assert(correct_st <= correct_ed);

    /* evaluate front end */
    {
      for (auto i = 0u; i < correct_st; i++) {
        if (cigar[i].op == 'I') {
          f_insertion += cigar[i].size;
        } else if (cigar[i].op == 'D') {
          f_deletion += cigar[i].size;
        } else if (cigar[i].op == 'X') {
          f_mismatch += cigar[i].size;
        } else if (cigar[i].op == '=') {
          f_match += cigar[i].size;
        }
      }
    }
    /* evaluate back end */
    {
      for (auto i = correct_ed + 1; i < cigar.size(); i++) {
        if (cigar[i].op == 'I') {
          b_insertion += cigar[i].size;
        } else if (cigar[i].op == 'D') {
          b_deletion += cigar[i].size;
        } else if (cigar[i].op == 'X') {
          b_mismatch += cigar[i].size;
        } else if (cigar[i].op == '=') {
          b_match += cigar[i].size;
        }
      }
    }
    /* evaluate middle part */
    {
      auto consecutive_lens = std::vector<std::size_t>{};
      auto no_del_len = 0ul;
      for (auto i = correct_st; i <= correct_ed; i++) {
        if (cigar[i].op != 'D') {
          no_del_len += cigar[i].size;
        } else {
          consecutive_lens.push_back(no_del_len);
        }
      }
      consecutive_lens.push_back(no_del_len);
      if (!consecutive_lens.empty()) {
        avg_del_len =
            std::accumulate(consecutive_lens.begin(), consecutive_lens.end(),
                            0.0) /
            consecutive_lens.size();
      }
    }
    // auto s = std::string{};
    // for (auto i = correct_st; i <= correct_ed; i++) {
    //   s += std::to_string(cigar[i].size) + cigar[i].op;
    // }
    // spdlog::debug("corrected cigar: {}", s);
  }

  /**
   * @brief evaluation on variant that covered by read
   *
   * @tparam R range for VcfRecord
   * @param vcf std::range::subrange for VcfRecord, which contains variants like
   * snps or INDELs
   * @return None
   * ? is it work for indel or just snp
   */
  template <std::ranges::range R>
  auto variant_eval(R vcf, bio::MafRecord& maf) {

    // auto print_seq = [](std::string_view s, std::size_t pos) {
    //   const auto len = 25;
    //   auto st = std::max(0ul, pos - len);
    //   auto ed = std::min(s.size(), pos + len);
    //   spdlog::debug("{:6} : {}", st, s.substr(st, ed - st));
    // };

    auto offset = cal_offset(cigar);
    covered_variants += std::ranges::size(vcf);

    /**
     * get variant start and end related on haplotype coordinate system, and
     * change all coordinate to 0-based
     */
    for (const auto& v : vcf) {
      auto vinfo = parse_vcf_info(v);
      auto st = std::stoul(std::string(vinfo["sim_start"])) - 1 - maf.ref.start;
      auto ed = std::stoul(std::string(vinfo["sim_end"])) - 1 - maf.ref.start;
      st = transform_coordinate(offset, st);
      ed = transform_coordinate(offset, ed);
      if (st >= corrected_read.seq.size()) {
        v_out_of_range += 1;
        continue;
      }
      auto alt_in_cread = corrected_read.seq.substr(st, ed - st + 1);
      if (v.alt != alt_in_cread) {
        v_mismatch += 1;
      }
    }
  }

  // auto print() {
  //   std::cout << std::fixed << std::setprecision(4) << std::endl;
  //   std::cout << "correction result: " << std::endl;
  //   std::cout << "Use elector please" << std::endl << std::endl;
  //   // std::cout << "accuracy: " << mat.accuracy() << std::endl;
  //   // std::cout << "precision: " << mat.precision() << std::endl;
  //   // std::cout << "senstivity: " << mat.sensitivity() << std::endl;
  //   // std::cout << "specificity: " << mat.specificity() << std::endl <<
  //   // std::endl;

  //   std::cout << "variant correction result: " << std::endl;
  //   std::cout << "v_mismatch: " << v_mismatch << std::endl;
  //   std::cout << "v_out_of_range: " << v_out_of_range << std::endl;
  //   std::cout << "total_variants: " << total_covered_variants << std::endl;
  //   std::cout << "variant correct successful rate: "
  //             << 1 - (double) (v_mismatch + v_out_of_range) /
  //                        total_covered_variants
  //             << std::endl;
  // }

  std::string ref;
  bio::FastqRecord<false> raw_read;
  bio::FastaRecord<false> corrected_read;
  bio::Cigar cigar;

  /* whether the cigar string is a mess */
  bool has_corrected = true;

  /* variant metrics */
  ConfusionMatrix<> mat;
  std::size_t v_mismatch = 0u;
  std::size_t v_out_of_range = 0u;
  std::size_t covered_variants = 0u;

  /* front end correction metrics */
  std::size_t f_match = 0u;
  std::size_t f_insertion = 0u;
  std::size_t f_deletion = 0u;
  std::size_t f_mismatch = 0u;

  /* back end correction metrics */
  std::size_t b_match = 0u;
  std::size_t b_insertion = 0u;
  std::size_t b_deletion = 0u;
  std::size_t b_mismatch = 0u;

  /* correction metrics except front end and back end of the read */
  constexpr static auto MAX_STATE_LEN = 10u;
  std::array<std::size_t, MAX_STATE_LEN> e_insertion = {};
  std::array<std::size_t, MAX_STATE_LEN> e_deletion = {};
  std::array<std::size_t, MAX_STATE_LEN> e_mismatch = {};
  double avg_del_len = 0.0;
};

template<class Head, class... Tail> 
auto add_one_row(rapidcsv::Document& doc, int row, int col, const Head& h, const Tail&... tails) {
  if constexpr (std::same_as<Head, std::string>) {
    doc.SetCell<Head>(col, row, h);
  } else {
    doc.SetCell<std::string>(col, row, std::to_string(h));
  }
  if constexpr (sizeof...(tails) > 0) {
    add_one_row(doc, row, col + 1, tails...);
  }
}

auto save_result(
    const std::vector<CorrectedResult>& results, const fs::path& path) {
  if (path.extension() != ".csv") {
    spdlog::error("Output file must be csv format, skip saving result");
    return;
  }
  rapidcsv::Document doc(std::string{}, rapidcsv::LabelParams(0, -1),
                         rapidcsv::SeparatorParams(',', true));  
  
  auto columnName = std::vector<std::string>{
    "name", "has_corrected", 
    "raw_len", "corrected_len", "perfect_len", 
    "cigar",
    "v_mismatch", "v_out_of_range", "covered_variants",
    "f_match", "f_insertion", "f_deletion", "f_mismatch",
    "b_match", "b_insertion", "b_deletion", "b_mismatch",
    "e_insertion", "e_deletion", "e_mismatch", "avg_del_len"
  };
  for (auto i = 0u; i < columnName.size(); i++) {
    doc.SetColumnName(i, columnName[i]);
  }
  for (auto& res : results) {
    add_one_row(doc, doc.GetRowCount(), 0, 
      res.raw_read.name, res.has_corrected, 
      res.raw_read.seq.size(), res.corrected_read.seq.size(), res.ref.size(),
      res.has_corrected ? std::string(res.cigar) : "",
      res.v_mismatch, res.v_out_of_range, res.covered_variants,
      res.f_match, res.f_insertion, res.f_deletion, res.f_mismatch,
      res.b_match, res.b_insertion, res.b_deletion, res.b_mismatch,
      std::accumulate(res.e_insertion.begin(), res.e_insertion.end(), 0u),
      std::accumulate(res.e_deletion.begin(), res.e_deletion.end(), 0u),
      std::accumulate(res.e_mismatch.begin(), res.e_mismatch.end(), 0u),
      res.avg_del_len);

      // TODO: show INDELs and mismatch in middle part by using better metrics
  }
  spdlog::info("write result to {}", path.string());
  doc.Save(path.string());
}