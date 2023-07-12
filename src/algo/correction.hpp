#pragma once

#include <algorithm>
#include <execution>
#include <map>
#include <ranges>
#include <string>
#include <thread>
#include <vector>

#include <biovoltron/file_io/all.hpp>
#include <spdlog/spdlog.h>

/* ugly solution, fix it later */
#include "../utility/all.hpp"

class BaseReadCorrector {

private:
  /**
   * @brief
   *
   * @param reads
   * @return auto
   */
  auto read_preprocess(const std::vector<bio::FastaRecord<false>>& reads) {
    const auto min_read_length = 1000ull;
    auto filtered_reads = std::vector<bio::FastaRecord<false>>{};
    /* there's no std::ranges::move_if() */
    for (auto& r : reads) {
      if (r.seq.size() >= min_read_length) {
        filtered_reads.emplace_back(std::move(r));
      }
    }
    std::sort(std::execution::par, filtered_reads.begin(), filtered_reads.end(),
              [&](auto& a, auto& b) { return a.name < b.name; });
    // std::ranges::sort(filtered_reads, {}, &bio::FastaRecord<false>::name)
    return filtered_reads;
  }

  /**
   * @brief
   * @param raw_reads
   * @param overlaps
   * @return
   */
  auto overlap_preprocess(const std::vector<bio::FastaRecord<false>>& raw_reads,
                          std::vector<bio::PafRecord>& overlaps) {
    assert(
        std::ranges::is_sorted(raw_reads, {}, &bio::FastaRecord<false>::name) &&
        "raw_reads must be sorted by name");
    auto valid_overlap = [&](const bio::PafRecord& overlap) {
      auto valid_read = [&](const std::string_view name) {
        return std::ranges::binary_search(raw_reads, name, {},
                                          &bio::FastaRecord<false>::name);
      };
      bool valid = valid_read(overlap.qname) && valid_read(overlap.tname);
      /* filter gogogo */
      return valid;
    };
    const auto min_overlap_length = 1000ull;
    const auto min_overlap_identity = 0.8;
    /* there's no std::ranges::move_if() */
    auto filtered_overlaps = std::vector<bio::PafRecord>{};
    for (auto& overlap : overlaps) {
      if (valid_overlap(overlap)) {
        filtered_overlaps.emplace_back(std::move(overlap));
      }
    }
    overlaps_size = filtered_overlaps.size();
    std::sort(std::execution::par, filtered_overlaps.begin(),
              filtered_overlaps.end());
    return filtered_overlaps;
  }

public:
  BaseReadCorrector(std::vector<bio::FastaRecord<false>>&& raw_reads,
                    std::vector<bio::PafRecord>&& overlaps,
                    const std::string& platform,
                    const int thread_num = std::thread::hardware_concurrency())
      : raw_reads(raw_reads), overlaps(overlaps), platform(platform) {
    
    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();
    this->raw_reads = read_preprocess(raw_reads);
    raw_read_size = this->raw_reads.size();
    this->overlaps = overlap_preprocess(this->raw_reads, overlaps);
    /* transform read name to id */
    for (std::size_t i = 0; i < this->raw_reads.size(); ++i) {
      name2id[this->raw_reads[i].name] = i;
    }
  }

  auto correct() { spdlog::info("Base corrector don't correct read :("); }

  auto print_info() {
    spdlog::info("Total raw reads: {}", unfiltered_raw_read_size);
    spdlog::info("Total overlaps: {}", unfiltered_overlap_size);
    spdlog::info("Platform: {}", platform);
    spdlog::info("Filtered raw reads: {}", raw_read_size);
    spdlog::info("Filtered overlaps: {}", overlaps_size);
  }
protected:

  std::size_t unfiltered_raw_read_size;
  std::size_t unfiltered_overlap_size;

  std::vector<bio::FastaRecord<false>> raw_reads;
  std::vector<bio::PafRecord> overlaps;
  std::string platform;
  std::size_t raw_read_size;
  std::size_t overlaps_size;

  /* need id to name? */
  std::map<std::string_view, std::size_t> name2id;
};

class FragmentedReadCorrector : public BaseReadCorrector {
public:
  using BaseReadCorrector::BaseReadCorrector;

  auto correct() {
    /* O(n) get overlap on all reads, but we can do O(n*log(n)) for convinent */

    auto get_overlap_range = [&](std::string_view read_name) {
      static auto idxR = 0u;
      auto st = idxR;
      while (idxR < overlaps.size() && overlaps[idxR].qname == read_name) {
        ++idxR;
      }
      // assert(idxR > st && "No overlap on this read");
      return std::ranges::subrange(overlaps.begin() + st,
                                   overlaps.begin() + idxR);
    };

    auto mn = std::numeric_limits<std::size_t>::max();
    auto mx = std::numeric_limits<std::size_t>::min();
    auto cnt = 0u;

    // /* try using spoa */

    // auto alignment_engine =
    //     spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);
    // auto graph = spoa::Graph{};

    // for (auto& seq : s) {
    //   auto alignment = alignment_engine->Align(seq, graph);
    //   graph.AddAlignment(alignment, seq);
    // }

    // auto consensus = graph.GenerateConsensus();
    // spdlog::info("consensus = {}", consensus);

    // auto msa = graph.GenerateMultipleSequenceAlignment();
    // for (auto& it : msa) {
    //   spdlog::info("msa seq = {}", it);
    // }

    for (auto& raw_read : raw_reads | std::views::take(2)) {
      spdlog::info("read name: {}", raw_read.name);
      auto range = get_overlap_range(raw_read.name);
      // TODO: extend overlap

      for (const auto& overlap : range) {
        /* split overlap into window */
        /* what inside a window */
        // 1. the subsequence of query read
        // 2. the overlap part of the subsequence from target reads
        
        auto ss = std::stringstream{};
        ss << overlap;
        spdlog::info("overlap: {}", ss.str());
        // spdlog::info("t.start = {}, t.end = {}, q.start = {}, q.end = {}",
        //              overlap.tstart, overlap.tend, overlap.qstart,
        //              overlap.qend);
      }
      auto sz = std::ranges::size(range);
      if (sz == 0) {
        cnt += 1;
        continue;
      }
      mn = std::min(mn, sz);
      mx = std::max(mx, sz);
      
      // spdlog::info("overlap size: {}", std::ranges::ssize(range));
      // for (auto& overlap : get_overlap_range(raw_read.name)) {
      //   spdlog::info("{} <-> {}", overlap.tname, overlap.qname);
      // }
    }
    spdlog::info("no overlap read: {}", cnt);
    spdlog::info("min overlap size: {}", mn);
    spdlog::info("max overlap size: {}", mx);

    auto corrected_read = std::vector<bio::FastaRecord<false>>{};
    return corrected_read;
  }

private:
  std::vector<std::vector<std::string>> read_fragments;
  /* maybe segment tree or faster data structure to calcualte coverage */
};