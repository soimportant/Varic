#pragma once

#include <concepts>
#include <execution>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include <biovoltron/file_io/fasta.hpp>

#include "thesis/format/paf.hpp"

namespace bio = biovoltron;

template <class R>
  requires std::derived_from<R, bio::FastaRecord<false>>
class BaseReadCorrector {
 private:
  /**
   * @brief
   *
   * @param reads
   * @return auto
   */
  auto read_preprocess() {
    spdlog::info("Read preprocessing...");
    const auto min_read_length = 500ull;
    auto filtered_reads = std::vector<R>{};
    /* there's no std::ranges::move_if() */
    for (auto& r : this->raw_reads) {
      if (r.seq.size() >= min_read_length) {
        filtered_reads.emplace_back(std::move(r));
      }
    }
    std::sort(std::execution::par, filtered_reads.begin(), filtered_reads.end(),
              [&](auto& a, auto& b) { return a.name < b.name; });
    raw_reads_size = filtered_reads.size();
    reads.resize(raw_reads_size);
    for (auto i : std::views::iota(0ull, raw_reads_size)) {
      reads[i] =
    }


    std::swap(this->raw_reads, filtered_reads);

    /* transform read name to id */
    for (const auto& raw_read : this->raw_reads) {
      name2id[raw_read.name] = name2id.size();
    }
    reads.resize(raw_reads_size);
  }

  /**
   * @brief
   * @param raw_reads
   * @param overlaps
   * @return
   */
  auto overlap_preprocess(std::vector<bio::PafRecord>& overlaps) {
    assert(std::ranges::is_sorted(raw_reads, {}, &R::name) &&
           "raw_reads must be sorted by name");
    spdlog::info("Overlap preprocessing...");

    const auto min_overlap_length = 500ull;
    // In original .paf file of overlaps between raw_reads, the match base
    // devided by alignment length is very low, figure out how it be computed
    // and why it is so low.
    // const auto min_overlap_identity = 0.7l;

    auto valid_overlap = [&](const bio::PafRecord& overlap) {
      auto valid_read = [&](const std::string_view name) {
        return std::ranges::binary_search(raw_reads, name, {}, &R::name);
      };
      if (!valid_read(overlap.q_name) || !valid_read(overlap.t_name)) {
        return false;
      }
      if (overlap.q_name == overlap.t_name) {
        return false;
      }
      if (overlap.aln_len < min_overlap_length) {
        return false;
      }
      {
        // filter conditions used by racon
        // 1. the length of alignment block from query read and target read
        // should not be too different
        // 2. remove internal match
        //  - A match is internal match if
        //      (overhang length > match length * internal threshold)
        //  - internal threshold = 0.8 for racon by using tools fpa
        //    -
        //    https://github.com/natir/fpa/blob/master/src/filter/internalmatch.rs

        auto q_overlap_len = overlap.q_end - overlap.q_start + 1;
        auto t_overlap_len = overlap.t_end - overlap.t_start + 1;
        auto possible_match_ratio =
            std::min(q_overlap_len, t_overlap_len) /
            (double) std::max(q_overlap_len, t_overlap_len);
        if (possible_match_ratio < 0.7) {
          return false;
        }

        const double internal_threshold = 0.8;
        auto q_head_len = overlap.q_start;
        auto q_tail_len = overlap.q_len - overlap.q_end;
        auto t_head_len = overlap.t_start;
        auto t_tail_len = overlap.t_len - overlap.t_end;
        auto overhang_len =
            overlap.strand == '+'
                ? std::min(q_head_len,
                           t_head_len + std::min(q_tail_len, t_tail_len))
                : std::min(q_head_len,
                           t_tail_len + std::min(q_tail_len, t_head_len));
        auto alignment_block_len = std::max(q_overlap_len, t_overlap_len);
        if (overhang_len >
            std::min(1000.0, alignment_block_len * internal_threshold)) {
          return false;
        }
      }

      // ? filter by identity
      // if (overlap.aln_len * min_overlap_identity > overlap.match) {
      //   return false;
      // }
      return true;
    };

    auto filtered_overlaps = std::vector<bio::PafRecord>{};

    /* there's no std::ranges::move_if(), so sad */
    // TODO: parallel do here, but need to care about data race when creating
    // TODO: reverse complement sequence and quality
    // TODO: may use another array for marking which read need to create rc
    for (auto& overlap : overlaps) {
      if (valid_overlap(overlap)) {
        /* add reverse complement sequence and reverse quality if needed */
        if (overlap.strand == '-') {
          auto tid = name2id[overlap.t_name];
          if (raw_reads_rc_seq[tid].empty()) {
            raw_reads_rc_seq[tid] = bio::Codec::rev_comp(raw_reads[tid].seq);
            if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
              raw_reads_rev_qual[tid] = raw_reads[tid].qual;
              std::ranges::reverse(raw_reads_rev_qual[tid]);
            }
          }
        }
        filtered_overlaps.emplace_back(std::move(overlap));
      }
    }
    overlaps_size = filtered_overlaps.size();
    std::sort(std::execution::par, filtered_overlaps.begin(),
              filtered_overlaps.end());
    // auto tmp_file = "/mnt/ec/ness/yolkee/thesis/tests/tmp/filtered_overlap.paf";
    // auto fout = std::ofstream(tmp_file);
    // for (const auto& overlap : filtered_overlaps) {
    //   fout << overlap << std::endl;
    // }
    std::swap(this->overlaps, filtered_overlaps);
  }

 protected:
  BaseReadCorrector(std::vector<R>&& raw_reads,
                    std::vector<bio::PafRecord>&& overlaps,
                    const std::string& platform,
                    const int thread_num = std::thread::hardware_concurrency(),
                    bool debug = false)
      : raw_reads(raw_reads),
        overlaps(overlaps),
        platform(platform),
        threads(thread_num),
        debug(debug) {
    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();

    read_preprocess();
    overlap_preprocess();
    print_info();
  }

  auto correct() -> void {
    spdlog::warn("Base corrector don't correct read :(");
  }

  auto print_info() -> void {
    spdlog::info("Total raw reads: {}", unfiltered_raw_read_size);
    spdlog::info("Total overlaps: {}", unfiltered_overlap_size);
    spdlog::info("Platform: {}", platform);
    spdlog::info("Filtered raw reads: {}", raw_reads_size);
    spdlog::info("Filtered overlaps: {}", overlaps_size);
    spdlog::info("Debug mode is {}", debug ? "on" : "off");
  }

  /* how much reads in original data */
  std::size_t unfiltered_raw_read_size;
  /* how much overlap records in original data */
  std::size_t unfiltered_overlap_size;

  /* raw TGS long reads */
  std::vector<R> raw_reads;
  std::vector<std::string> raw_reads_rc_seq;
  std::vector<std::string> raw_reads_rev_qual;

  /* overlap information between `raw_reads` */
  std::vector<bio::PafRecord> overlaps;

  std::size_t raw_reads_size;
  std::size_t overlaps_size;

  /* sequencing platform of raw_reads */
  // TODO: should be enum class
  std::string platform;

  /* need id to name? */
  std::map<std::string, std::size_t> name2id;

  /* threads */
  std::size_t threads;

  /* debug flag */
  bool debug = false;
};