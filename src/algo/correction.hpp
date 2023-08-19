#pragma once

#include <algorithm>
#include <concepts>
#include <execution>
#include <map>
#include <ranges>
#include <string>
#include <thread>
#include <vector>
#include <random>

#include <biovoltron/file_io/all.hpp>
#include <spdlog/spdlog.h>

// ! ugly solution, fix it later
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
    spdlog::info("Read proprecessing...");
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
    raw_reads_rc_seq.resize(filtered_reads.size());
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
    spdlog::info("Overlap proprecessing...");

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
    auto filtered_overlaps = std::vector<bio::PafRecord>{};

    /* there's no std::ranges::move_if(), so sad */
    for (auto& overlap : overlaps) {
      if (valid_overlap(overlap)) {
        /* add reverse complement if we need */
        if (overlap.strand == '-') {
          auto tid = name2id[overlap.tname];
          if (raw_reads_rc_seq[tid].empty()) {
            raw_reads_rc_seq[tid] = bio::Codec::rev_comp(raw_reads[tid].seq);
          }
        }
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
                    const int thread_num = std::thread::hardware_concurrency(),
                    bool debug = false)
      : raw_reads(raw_reads), overlaps(overlaps), platform(platform) {

    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();
    this->raw_reads = read_preprocess(raw_reads);
    raw_read_size = this->raw_reads.size();
    /* transform read name to id */
    for (std::size_t i = 0; i < this->raw_reads.size(); ++i) {
      name2id[this->raw_reads[i].name] = i;
    }
    this->overlaps = overlap_preprocess(this->raw_reads, overlaps);
    overlaps_size = this->overlaps.size();
    threads = thread_num;
    this->debug = debug;
    print_info();
  }

  auto correct() { spdlog::warn("Base corrector don't correct read :("); }

  auto print_info() -> void {
    spdlog::info("Total raw reads: {}", unfiltered_raw_read_size);
    spdlog::info("Total overlaps: {}", unfiltered_overlap_size);
    spdlog::info("Platform: {}", platform);
    spdlog::info("Filtered raw reads: {}", raw_read_size);
    spdlog::info("Filtered overlaps: {}", overlaps_size);
    spdlog::info("Debug mode is {}", debug ? "on" : "off");
  }

protected:
  /* how much reads in original data */
  std::size_t unfiltered_raw_read_size;
  /* how much overlap records in original data */
  std::size_t unfiltered_overlap_size;

  /* raw TGS long reads */
  std::vector<bio::FastaRecord<false>> raw_reads;
  std::vector<std::string> raw_reads_rc_seq;

  /* overlap information between `raw_reads` */
  std::vector<bio::PafRecord> overlaps;

  std::size_t raw_read_size;
  std::size_t overlaps_size;

  /* sequencing platform of raw_reads */
  // TODO: should be enum class
  std::string platform;

  /* need id to name? */
  std::map<std::string, std::size_t> name2id;

  /* threads */
  std::size_t threads;

  /* debug */
  bool debug = false;
};

class FragmentedReadCorrector : public BaseReadCorrector {
public:
  using BaseReadCorrector::BaseReadCorrector;

  class Window {
  public:
    // basic length of a window
    static const auto window_len = 500ul;
    // overlap length of adjanency window. Therefore, typeical length of a
    // single window should be window_len + 2 * window_overlap_len
    static const auto window_overlap_len = 50ul;
    // extend length of overlap region
    static const auto overlap_extend_len = 25ul;

    /* backbone read information */
    std::size_t idxL, idxR;
    std::string backbone_seq;

    /* overlap information of backbone read */
    std::vector<std::size_t> overlap_read_idx;
    std::vector<std::pair<std::size_t, std::size_t>> overlap_boundaries;
    std::vector<std::string> overlap_seqs;
  };

  class Read : public bio::FastaRecord<false> {
  public:
    // coverage -> segment tree like data structure

    // collect correct seqs from other reads
    // may further assemble to correct version of this read
    std::size_t idx;
    std::string rc_seq;
    std::vector<std::string> corrected_fragments;

    std::vector<Window> windows;
  };

  /**
   * @brief Get the overlap range object
   *
   * @param read_name
   * @return std::range::range
   */
  std::ranges::range auto get_overlap_range(const std::string_view read_name) {
    static auto idxR = 0u;
    auto st = idxR;
    while (idxR < overlaps.size() && overlaps[idxR].qname == read_name) {
      ++idxR;
    }
    // assert(idxR > st && "No overlap on this read");
    return std::ranges::subrange(overlaps.begin() + st,
                                 overlaps.begin() + idxR);
  }

  /**
   * @brief Generate windows from the overlap between raw_read and other reads.
   *
   * @tparam R constrain by std::ranges::range
   * @param raw_read
   * @param overlap_range
   * @return std::vector<Window>
   */
  auto make_windows(const Read& raw_read) {

    auto windows = std::vector<Window>{};

    /* get all overlap of this read */
    for (const auto& overlap : get_overlap_range(raw_read.name)) {
      std::stringstream ss;
      ss << overlap;

      /* target boundaries information */
      auto boundary = std::vector<
          std::pair<std::size_t, std::pair<std::size_t, std::size_t>>>{};
      auto tid = name2id[overlap.tname];

      /* extend overlap boundaries at query read and target read */
      auto qstart = overlap.qstart > Window::overlap_extend_len
                        ? overlap.qstart - Window::overlap_extend_len
                        : 0ul;
      auto qend =
          std::min(overlap.qlen - 1, overlap.qend + Window::overlap_extend_len);
      auto qlen = qend - qstart + 1;

      auto tstart = overlap.tstart > Window::overlap_extend_len
                        ? overlap.tstart - Window::overlap_extend_len
                        : 0ul;
      auto tend =
          std::min(overlap.tlen - 1, overlap.tend + Window::overlap_extend_len);
      auto tlen = tend - tstart + 1;

      /* set initial boundary */
      for (auto widx = qstart / Window::window_len;
           widx * Window::window_len < qend; ++widx) {
        auto qidxL = std::max(qstart, widx * Window::window_len);
        auto qidxR =
            std::min((widx + 1) * Window::window_len - 1, overlap.qlen - 1);
        auto tidxL = std::min(tstart + (qidxL - qstart), overlap.tlen - 1);
        auto tidxR = tstart + (qidxR - qstart);
        boundary.emplace_back(widx, std::make_pair(tidxL, tidxR));
      }

      /* extend boundaries for making adjanency window have overlap */
      std::ranges::for_each(boundary, [&](auto& b) {
        auto& [_, tidx] = b;
        auto& [tidxL, tidxR] = tidx;
        tidxL = tidxL > Window::window_overlap_len
                    ? tidxL - Window::window_overlap_len
                    : 0;
        assert(tidxL < reads[tid].seq.size() && "tidxL out of range");
        tidxR = std::min(tidxR + Window::window_overlap_len, overlap.tlen - 1);
      });

      /* adjust boundaries by length difference */
      auto len_diff = (std::int64_t) (qlen - tlen);
      auto offset = len_diff / (std::int64_t) boundary.size();
      for (auto i = 1; i < (int) boundary.size(); i++) {
        boundary[i - 1].second.second += offset * i;
      }
      boundary.back().second.second =
          std::min(boundary.back().second.second + len_diff, tlen - 1);

      /* push sequence into window */
      for (const auto& [widx, tidx] : boundary) {
        if (widx >= windows.size()) {
          windows.resize(widx + 1);
        }
        auto tseq = std::string_view{};
        tseq = (overlap.strand == '+') ? reads[tid].seq
                                       : reads[tid].rc_seq;
        auto [tidxL, tidxR] = tidx;

        // TODO: set a filter here, for the sequence which length is less than
        // TODO: some value here, drop it.

        windows[widx].overlap_read_idx.push_back(tid);
        windows[widx].overlap_boundaries.emplace_back(tidxL, tidxR);
        windows[widx].overlap_seqs.emplace_back(
            tseq.substr(tidxL, tidxR - tidxL + 1));
      }
    }

    // ? if the overlap strand == '-', should we keep this information for
    // ? assemble the correcetd read, since the strand == '-', that means the
    // ? consensus produce by this read, would be the reverse complement of 
    // ? the target read, then we need to use this instead.

    /* set backbone sequence */
    for (auto widx : std::views::iota(0u, windows.size())) {
      auto idxL = widx * Window::window_len;
      auto idxR = (widx + 1) * Window::window_len - 1;
      if (idxL > Window::window_overlap_len) {
        idxL -= Window::window_overlap_len;
      }
      idxR =
          std::min(idxR + Window::window_overlap_len, raw_read.seq.size() - 1);
      windows[widx].idxL = idxL;
      windows[widx].idxR = idxR;
      windows[widx].backbone_seq = raw_read.seq.substr(idxL, idxR - idxL + 1);
    }

    return windows;
  }

  auto make_reads() {
    reads.resize(raw_read_size);
    for (auto i : std::views::iota(0ul, raw_read_size)) {
      reads[i].name = std::move(raw_reads[i].name);
      reads[i].seq = std::move(raw_reads[i].seq);
      reads[i].rc_seq = std::move(raw_reads_rc_seq[i]);
      reads[i].idx = name2id[reads[i].name];
      /* initial read */
    }
    // std::mt19937 rng(42);
    // std::ranges::shuffle(reads, rng);
  }

  auto correct() {
    make_reads();
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

    auto takes = debug ? 1 : std::numeric_limits<std::size_t>::max();
    
    for (auto& raw_read : reads | std::views::take(takes)) {
      auto windows = make_windows(raw_read);
      spdlog::info("read name = {}, seq size = {}", raw_read.name,
                   raw_read.seq.size());
      spdlog::info("windows.size() = {}", windows.size());

      auto seq_cnt = 0;
      for (const auto& w : windows | std::views::take(takes)) {
        seq_cnt += w.overlap_seqs.size();
        
        // 1. build de burjin graph
        // 2. msa then use spoa graph
        // try to find some pruning policy for doing this.
        // + prune without realignment
        // 
        {
          auto alignment_engine =
              spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5,
                                            -3);
          auto graph = spoa::Graph{};
          for (const auto &seq : w.overlap_seqs) {
            auto alignment = alignment_engine->Align(seq, graph);
            graph.AddAlignment(alignment, seq);
          }

        }

        
        /**
         * 
         * 
         */

        /* do msa here */
        // auto alignment_engine =
        //     spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5,
        //     -3);
        // auto graph = spoa::Graph{};
        // auto alignment = alignment_engine->Align(w.overlap_seqs[0], graph);
        // for (auto& seq : w.overlap_seqs | std::views::drop(1)) {
        //   auto s = std::string(seq);
        //   auto alignment = alignment_engine->Align(s, graph);
        //   graph.AddAlignment(alignment, s);
          
        // }
        
        // auto consensus = graph.GenerateConsensus();
        // spdlog::debug("get consensus");
        // spdlog::debug("idxL = {}, idxR = {}", w.idxL, w.idxR);
        // spdlog::debug("size = {}", consensus.size());
        // spdlog::debug("consensus = {}", consensus);
        // auto msa = graph.GenerateMultipleSequenceAlignment();
        // for (auto& it : msa) {
        //   spdlog::info("{:4} msa seq = {}", it.size(), it.substr(0, 100));
        // }
      }
      spdlog::debug("read name = {}, avg seq count = {}", raw_read.name, (double) seq_cnt / windows.size());
    }

    // for (auto& raw_read : raw_reads | std::views::take(1)) {
    //   spdlog::info("read name: {}", raw_read.name);
    //   // auto windows = std::vector<Window>();

    //   /* TODO: assign backbone sequence to each window */

    //   /* get all overlap of this read */
    //   auto range = get_overlap_range(raw_read.name);
    //   for (const auto& overlap : range) {
    //     /* adjust query sequence range */
    //     auto qstart = overlap.qstart > overlap_extend_len
    //                       ? overlap.qstart - overlap_extend_len
    //                       : 0ul;
    //     auto qend =
    //         std::min(overlap.qlen - 1, overlap.qend + overlap_extend_len);
    //     /* adjust target sequence range */
    //     auto tstart = overlap.tstart > overlap_extend_len
    //                       ? overlap.tstart - overlap_extend_len
    //                       : 0ul;
    //     auto tend =
    //         std::min(overlap.tlen - 1, overlap.tend + overlap_extend_len);

    //     auto widx = qstart / window_len;
    //     auto tid = name2id[overlap.tname];

    //     spdlog::debug("tlen = {}, qlen = {}, diff = {}", tend - tstart,
    //                   qend - qstart,
    //                   (std::int64_t)(tend - tstart - (qend - qstart)));
    //     for (auto widx = qstart / window_len;
    //          widx * window_len + window_overlap_len < qend; widx++) {
    //       auto qidxL = std::max(
    //           qstart, widx * window_len - (widx == 0 ? 0 :
    //           window_overlap_len));
    //       auto qidxR =
    //           std::min(qend, (widx + 1) * window_len + window_overlap_len);
    //       // The sequence length of query part and target part may not
    //       // equal, maybe we should dynamic adjust, or extend or shrink
    //       // some random length or some value according to the difference
    //       // between query part and target part.
    //       auto tidxL = std::max(tstart, overlap.tstart + qidxL - qstart);
    //       auto tidxR = std::min(tend, overlap.tstart + qidxR - qstart);

    //       // spdlog::debug("qidxL = {}, qidxR = {}, tidxL = {}, tidxR = {}",
    //       // qidxL,
    //       //               qidxR, tidxL, tidxR);

    //       auto tseq = std::string_view{};
    //       tseq = (overlap.strand == '+') ? raw_reads[tid].seq
    //                                      : raw_reads_rc_seq[tid];
    //       window[widx].emplace_back(tseq.substr(tidxL, tidxR - tidxL));
    //     }
    //   }

    //   /* for each overlap, assign target sequence to each window */
    //   // for (const auto& overlap : range) {
    //   //   auto widx = overlap.qstart / window_len;
    //   //   while ((widx + 1) * window_len + extend_len < overlap.qend) {
    //   //     auto tid = name2id[overlap.tname];
    //   //     auto qstart = std::max(overlap.qstart, widx * window_len -
    //   //     extend_len); auto qend = std::min(overlap.qend, (widx + 1) *
    //   //     window_len + extend_len);

    //   //     /* The sequence length of query part and target part may not
    //   equal
    //   //     */
    //   //     /* Maybe we should dynamic adjust, or extend or shrink random
    //   //     length */ auto tstart = overlap.tstart + qstart -
    //   overlap.qstart;
    //   //     auto tend = std::min(
    //   //       raw_reads[tid].seq.size(),
    //   //       overlap.tstart + qend - overlap.qstart
    //   //     );
    //   //     auto tseq = overlap.strand == '+' ? raw_reads[tid].seq :
    //   //                                         raw_reads_rc_seq[tid];
    //   //     tseq = tseq.substr(tstart, tend - tstart);
    //   //     window[name2id[overlap.qname]][widx].emplace_back(tseq);
    //   //     widx++;
    //   //   }
    //   // }

    //   // for (auto idxL = 0ull, idxR = 0ull; idxR < range.size() ; idxR++) {
    //   //   while (idxL < idxR && range[idxL].qlen < window_lbound) {
    //   //     idxL++;
    //   //   }
    //   //   std::vector<std::string_view> seqs;
    //   //   auto window_rbound = window_lbound + window_len;

    //   //   /* add backbone sequence into window */
    //   //   seqs.emplace_back(std::string_view(raw_read.seq).substr(
    //   //     window_lbound - extend_len,
    //   //     std::min(window_len + extend_len, raw_read.seq.size() -
    //   //     window_lbound)
    //   //   ));
    //   //   for (auto i = idxL; i < idxR; i++) {
    //   //     /* add overlap sequence into window */
    //   //     auto tid = name2id[range[i].tname];
    //   //     auto tseq = std::string_view{};
    //   //     auto tstart =
    //   //     if (range[i].strand == '+') {
    //   //       tseq =
    //   //       std::string_view(raw_reads[tid].seq).substr(range[i].tstart,
    //   //       range[i].tend);
    //   //     } else {

    //   //     }
    //   //   }
    //   //   window_lbound += window_len;
    //   // }

    //   // spdlog::info("overlap size: {}", std::ranges::ssize(range));
    //   // for (auto& overlap : get_overlap_range(raw_read.name)) {
    //   //   spdlog::info("{} <-> {}", overlap.tname, overlap.qname);
    //   // }
    // }

    spdlog::info("no overlap read: {}", cnt);
    spdlog::info("min overlap size: {}", mn);
    spdlog::info("max overlap size: {}", mx);

    auto corrected_read = std::vector<bio::FastaRecord<false>>{};
    return corrected_read;
  }
private:
  std::vector<Read> reads;
};