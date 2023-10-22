#pragma once

#include <edlib.h>

#include <concepts>
#include <future>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <ranges>
#include <spoa/spoa.hpp>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "thesis/algo/assemble/read_assembler.hpp"
#include "thesis/corrector/base_corrector.hpp"
#include "thesis/corrector/detail/read.hpp"
#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/corrector/detail/window.hpp"
#include "thesis/utility/threadpool/threadpool.hpp"

auto align(const std::string_view q, std::string_view t) {
  auto aln = edlibAlign(
      q.data(), q.size(), t.data(), t.size(),
      edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

  if (aln.status != EDLIB_STATUS_OK) {
    spdlog::error(
        "Error when align overlap part of query read and target read");
    exit(EXIT_FAILURE);
  }
  /* there's no *X* in cigar string when using `EDLIB_CIGAR_STANDARD */
  /* if you wanna contain mismatch info, please use `EDLIB_CIGAR_EXTENDED */
  auto* s = edlibAlignmentToCigar(aln.alignment, aln.alignmentLength,
                                  EDLIB_CIGAR_STANDARD);
  auto cigar = bio::Cigar(s);
  std::free(s);
  edlibFreeAlignResult(aln);
  return cigar;
}

/**
 * @brief following the implementation of `find_breakpoints_from_cigar` in racon
 * but do a little bit optimization
 *
 * @link https://github.com/isovic/racon/blob/master/src/overlap.cpp
 *
 */
template <std::ranges::range R>
auto find_breakpoints_from_cigar(const bio::Cigar& cigar, const R window_range,
                                 const std::size_t q_start,
                                 const std::size_t q_end,
                                 const std::size_t t_start) {
  auto window_cnt = std::ranges::size(window_range);

  // for (auto& w : window_range) {
  //   spdlog::debug("w.start = {}, w.end = {}", std::max(q_start, w.start),
  //   std::min(q_end, w.end));
  // }

  /* corresponding interval of target read for each window of query read */
  auto t_breakpoints =
      std::vector<std::pair<std::size_t, std::size_t>>(window_cnt);

  /* iterator for query read and target read */
  auto q_iter = q_start, t_iter = t_start;
  auto t_last_match = t_iter;

  /* store the window indexes that window is not begin at 'M' region */
  auto wait_for_first_match = std::queue<std::size_t>{};

  /* store the window indexes that not reach the end of window */
  auto wait_for_last_match = std::queue<std::size_t>{};

  auto w_idx = 0u;
  auto update = [&]() {
    while (!wait_for_last_match.empty()) {
      auto idx = wait_for_last_match.front();
      /* q_end < window_range[idx].end -> the last window */
      if (q_iter + 1 == std::min(q_end, window_range[idx].end)) {
        wait_for_last_match.pop();
        t_breakpoints[idx].second = t_last_match + 1;
      } else {
        break;
      }
    }
  };

  for (auto& [len, op] : cigar) {
    if (op == 'M') {
      while (!wait_for_first_match.empty()) {
        auto idx = wait_for_first_match.front();
        wait_for_first_match.pop();
        t_breakpoints[idx].first = t_iter;
        wait_for_last_match.push(idx);
        // spdlog::debug("[{}{}] push {}, q_iter = {}, t_iter = {}", len, op,
        // idx, q_iter,
        //               t_iter);
      }
      for (auto i = 0; i < len; i++, q_iter++, t_iter++) {
        /**
         * q_start > window_range[w_idx].start -> the first window
         *
         * while loop may run twice when two window start at same position,
         * which means the length of first window is less than
         * `Window::window_extend_len`
         */
        while (w_idx < window_cnt &&
               q_iter == std::max(q_start, window_range[w_idx].start)) {
          wait_for_last_match.push(w_idx);
          // spdlog::debug("push {}, q_iter = {}, t_iter = {}", w_idx, q_iter,
          // t_iter);
          t_breakpoints[w_idx].first = t_iter;
          w_idx++;
        }
        t_last_match = t_iter;
        update();
      }
    } else if (op == 'I') {
      for (auto i = 0; i < len; i++, q_iter++) {
        while (w_idx < window_cnt &&
               q_iter == std::max(q_start, window_range[w_idx].start)) {
          wait_for_first_match.push(w_idx);
          w_idx++;
        }
        update();
      }
    } else if (op == 'D') {
      t_iter += len;
    } else {
      spdlog::error("unknown cigar operation: {}", op);
      std::exit(EXIT_FAILURE);
    }
  }
  if (wait_for_last_match.size() != 0) {
    spdlog::debug("q_iter = {}, t_iter = {}", q_iter, t_iter);
    spdlog::debug("wait_for_last_match.size() = {}, front = {}",
                  wait_for_last_match.size(), wait_for_last_match.front());
    spdlog::debug(
        "wait_for_last_match.front().start = {}, "
        "wait_for_last_match.front().end = {}",
        window_range[wait_for_last_match.front()].start,
        window_range[wait_for_last_match.front()].end);
    for (auto& w : window_range) {
      spdlog::debug("window_st = {}, window_ed = {}",
                    std::max(q_start, w.start), std::min(q_end, w.end));
    }
    for (auto& [t_st, t_ed] : t_breakpoints) {
      spdlog::debug("t_st = {}, t_ed = {}", t_st, t_ed);
    }
  }
  assert(wait_for_last_match.empty());
  return t_breakpoints;
}

template <class R>
  requires std::derived_from<R, bio::FastaRecord<R::encoded>>
class FragmentedReadCorrector : public BaseReadCorrector<R> {
 public:
  /**
   * for the purpose of using `using BaseReadCorrector<R>::{}` here, please
   * refer to links below
   *
   * @link https://shorturl.at/hjqt3
   * @link https://shorturl.at/dnEM0
   *
   * The simplest explain is that we need to specify the data member type we
   * inherited from base class, otherwise, the compiler will not know which
   * type we are referring to.
   */
  using BaseReadCorrector<R>::BaseReadCorrector;
  using BaseReadCorrector<R>::overlaps;
  using BaseReadCorrector<R>::raw_reads;
  using BaseReadCorrector<R>::raw_reads_size;
  using BaseReadCorrector<R>::threads;
  using BaseReadCorrector<R>::name2id;
  using BaseReadCorrector<R>::debug;
  using BaseReadCorrector<R>::raw_reads_rc_seq;
  using BaseReadCorrector<R>::raw_reads_rev_qual;

  using Read = ReadWrapper<R>;

  /**
   * @brief Get the range of bio::PafRecord with the same read name.
   *
   * @param read_name name of the read
   * @return std::ranges::subrange
   */
  auto get_overlap_range(const std::string_view read_name) {
    auto st = std::ranges::lower_bound(overlaps, read_name, {},
                                       &bio::PafRecord::q_name);
    auto ed = std::ranges::upper_bound(overlaps, read_name, {},
                                       &bio::PafRecord::q_name);
    return std::ranges::subrange(st, ed);
  }

  /**
   * @brief
   *
   *
   * @tparam R constrained by std::ranges::range
   * @param raw_read
   * @param overlap_range
   * @return std::vector<Window>
   */
  auto make_windows(const Read& raw_read) {
    auto windows = std::vector<Window>{};

    /* initialize boundaries for each window */
    for (auto pos = 0ul; pos < raw_read.seq.size(); pos += Window::window_len) {
      auto w_idx = windows.size();
      auto window = Window{w_idx, pos, pos + Window::window_len};
      window.extend(raw_read.seq.size());
      std::tie(window.backbone.seq, window.backbone.qual) =
          raw_read.subview(true, window.start, window.end);
      // ? if the last window is too short, that means the right bound of the
      // ? second to last window is the end of read, then should we remove the
      // ? last window?
      windows.emplace_back(std::move(window));
    }

    auto get_window_range = [&](std::size_t start, std::size_t end) {
      auto st = std::ranges::lower_bound(windows, start, {}, &Window::start);
      auto ed = std::ranges::upper_bound(windows, end, {}, &Window::start);
      return std::ranges::subrange(st, ed);
    };

    for (const auto& overlap : get_overlap_range(raw_read.name)) {
      auto t_id = name2id[overlap.t_name];
      auto forward_strain = (overlap.strand == '+');

      auto [q_start, q_end] = std::tie(overlap.q_start, overlap.q_end);
      auto [t_start, t_end] = std::tie(overlap.t_start, overlap.t_end);

      // TODO: put this to PafRecord
      auto extend_overlap_range = [&]() {
        
        // done this way for avoiding unsigned integer underflow
        q_start =
            q_start > overlap_extend_len ? q_start - overlap_extend_len : 0ul;
        q_end = std::min(overlap.q_len, q_end + overlap_extend_len);
        t_start =
            t_start > overlap_extend_len ? t_start - overlap_extend_len : 0ul;
        t_end = std::min(overlap.t_len, t_end + overlap_extend_len);
      };
      extend_overlap_range();

      auto [q_seq_view, q_qual_view] = raw_read.subview(true, q_start, q_end);
      auto [t_seq_view, t_qual_view] =
          reads[t_id].subview(forward_strain, t_start, t_end);
      auto cigar = align(q_seq_view, t_seq_view);

      auto window_range = get_window_range(q_start, q_end);
      auto window_idx = window_range.begin() - windows.begin();
      auto breakpoints = find_breakpoints_from_cigar(cigar, window_range,
                                                     q_start, q_end, t_start);
      assert(std::ranges::size(window_range) == breakpoints.size() &&
             "size not equal between window_range and breakpoints");
      const auto sz = std::ranges::size(window_range);
      for (auto i = 0u; i < sz; i++, window_idx++) {
        auto& window = windows[window_idx];
        auto& [t_window_st, t_window_ed] = breakpoints[i];
        if (t_window_st >= t_end) {
          continue;
        }
        assert(t_window_ed >= t_window_st &&
               "t_window_ed < t_window_st, wrong boundary");
        auto [t_window_seq, t_window_qual] = reads[t_id].subview(
            forward_strain, t_window_st, t_window_ed);
        window.overlap_seqs.emplace_back(
            Sequence{.read_id = t_id,
                     .left_bound = t_window_st,
                     .right_bound = t_window_ed,
                     .seq = t_window_seq,
                     .qual = t_window_qual,
                     .forward_strain = forward_strain});
      }
    }

    // /* get all overlap of this read */
    // for (const auto& overlap : get_overlap_range(raw_read.name)) {
    //   /* boundaries information of target read -> (widx, (lbound, rbound)) */
    //   auto boundaries = std::vector<
    //       std::pair<std::size_t, std::pair<std::size_t, std::size_t>>>{};
    //   auto tid = name2id[overlap.tname];
    //   auto forward_strain = (overlap.strand == '+');

    //   auto qseq_view = std::string_view(raw_read.seq);
    //   auto qstart = overlap.qstart;
    //   auto qend = overlap.qend;

    //   auto tseq_view = std::string_view{};
    //   tseq_view = forward_strain ? reads[tid].seq : reads[tid].rc_seq;
    //   auto tstart =
    //       forward_strain ? overlap.tstart : overlap.tlen - overlap.tend;
    //   auto tend = forward_strain ? overlap.tend : overlap.tlen -
    //   overlap.tstart;

    //   /* extend overlap range at query read and target read */
    //   auto extend_overlap_range = [&]() {
    //     // done this way for avoiding unsigned integer underflow
    //     qstart = qstart > Window::overlap_extend_len
    //                  ? qstart - Window::overlap_extend_len
    //                  : 0ul;
    //     qend = std::min(overlap.qlen - 1, qend + Window::overlap_extend_len);
    //     tstart = tstart > Window::overlap_extend_len
    //                  ? tstart - Window::overlap_extend_len
    //                  : 0ul;
    //     tend = std::min(overlap.tlen - 1, tend + Window::overlap_extend_len);
    //   };

    //   /* initialize boundaries */
    //   auto init_boundaries = [&]() {
    //     for (auto widx = qstart / Window::window_len;
    //          widx * Window::window_len < qend; ++widx) {
    //       auto qidxL = std::max(qstart, widx * Window::window_len);
    //       auto qidxR = std::min((widx + 1) * Window::window_len - 1, qend);
    //       auto tidxL = std::min(tstart + (qidxL - qstart), tend);
    //       auto tidxR = std::min(tstart + (qidxR - qstart), tend);
    //       if (tidxL == tidxR) {
    //         continue;
    //       }
    //       boundaries.emplace_back(widx, std::make_pair(tidxL, tidxR));
    //     }
    //   };

    /**
     * adjust boundaries according to length difference between query
     * sequence and target sequence
     * ! this function may cause some boundary out of range
     * ! we may fix this by using the method from racon, that is,
     * ! align two overlaps, and determine w indow boundaries by cigar string
     */
    // auto adjust_boundaries = [&]() {
    //   auto len_diff = (std::int64_t)(tend - tstart + 1) -
    //                   (std::int64_t)(qend - qstart + 1);

    //   auto offset = len_diff / (std::int64_t) boundaries.size();
    //   for (auto i = 1u; i < boundaries.size(); i++) {
    //     auto& prev_tidxR = boundaries[i - 1].second.second;
    //     auto& curr_tidxL = boundaries[i].second.first;

    //     // do this way for avoiding unsigned integer underflow
    //     // I don't want debug this...
    //     if (offset >= 0) {
    //       prev_tidxR = std::min(prev_tidxR + offset * i, tend);
    //       curr_tidxL = std::min(curr_tidxL + offset * i, tend);
    //     } else {
    //       auto n_offset = -offset;
    //       prev_tidxR =
    //           (prev_tidxR > n_offset * i) ? prev_tidxR - n_offset * i :
    //           0ul;
    //       curr_tidxL =
    //           (curr_tidxL > n_offset * i) ? curr_tidxL - n_offset * i :
    //           0ul;
    //     }
    //   }
    //   auto& [last_tidxL, last_tidxR] = boundaries.back().second;
    //   last_tidxR = std::min(last_tidxR + len_diff, tend);
    //   if (last_tidxR <= last_tidxL) {
    //     boundaries.pop_back();
    //   }
    // };

    // /* extend boundaries for making adjanency window have overlap */
    // auto extend_boundaries = [&]() {
    //   std::ranges::for_each(boundaries, [&](auto& boundary) {
    //     auto& [_, tidx] = boundary;
    //     auto& [tidxL, tidxR] = tidx;
    //     /* cannot use std::min due to unsigned integer subtraction */
    //     tidxL = tidxL > Window::window_overlap_len
    //                 ? tidxL - Window::window_overlap_len
    //                 : 0ul;
    //     tidxR = std::min(tidxR + Window::window_overlap_len, tend);
    //   });
    // };

    // auto print = [&]() {
    //   spdlog::debug("qstart = {}, qend = {}, tstart = {}, tend = {}",
    //   qstart,
    //                 qend, tstart, tend);
    //   auto len_diff = (std::int64_t)(tend - tstart + 1) -
    //                   (std::int64_t)(qend - qstart + 1);
    //   auto offset = len_diff / (std::int64_t) boundaries.size();
    //   spdlog::debug("len_diff = {}, offset = {}", len_diff, offset);
    //   for (auto& [widx, tidx] : boundaries) {
    //     spdlog::debug("widx = {}, tidxL = {}, tidxR = {}", widx,
    //     tidx.first,
    //                   tidx.second);
    //   }
    // };

    // /* push sequence into window */
    // auto push_sequence = [&]() {
    //   for (const auto& [widx, tidx] : boundaries) {
    //     while (widx >= windows.size()) {
    //       windows.emplace_back(Window{});
    //     }
    //     auto forward_strain = (overlap.strand == '+');
    //     auto tseq = std::string_view{};
    //     auto tqual = std::optional<std::string_view>{};
    //     auto [tidxL, tidxR] = tidx;

    //     tseq = forward_strain ? reads[tid].seq : reads[tid].rc_seq;
    //     if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
    //       tqual = forward_strain ? reads[tid].qual : reads[tid].rc_qual;
    //     }

    //     if (tidxL > tidxR) {
    //       print();
    //     }
    //     assert(tidxL <= tidxR && "wrong boundary");
    //     assert(tidxL < tend && "tidxL out of range");

    //     auto subseq = tseq.substr(tidxL, tidxR - tidxL + 1);
    //     auto subqual = std::optional<std::string_view>{};
    //     if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
    //       subqual = tqual.value().substr(tidxL, tidxR - tidxL + 1);
    //     }
    //     assert(subseq.size() < Window::window_len * 2 &&
    //            "overlap sequence too long");
    //     if (subseq.size() < Window::overlap_min_len) {
    //       continue;
    //     }
    //     // TODO: Directly push sequence into graph
    //     windows[widx].overlaps.emplace_back(
    //         Overlap{.tid = tid,
    //                 .tidxL = tidxL,
    //                 .tidxR = tidxR,
    //                 .seq = subseq,
    //                 .qual = subqual,
    //                 .strain = forward_strain});
    //   }
    // };

    // extend_overlap_range();
    // init_boundaries();
    // // print();
    // adjust_boundaries();
    // // print();
    // extend_boundaries();
    // // print();
    // push_sequence();
    // }

    // ? if the overlap strand == '-', should we keep this information for
    // ? assemble the correcetd read, since the strand == '-', that means the
    // ? consensus produce by this read, would be the reverse complement of
    // ? the target read, then we need to use this instead.

    /* set backbone sequence */
    // for (auto widx : std::views::iota(0u, windows.size())) {
    //   auto idxL = widx * Window::window_len;
    //   auto idxR = (widx + 1) * Window::window_len - 1;
    //   if (idxL > Window::window_overlap_len) {
    //     idxL -= Window::window_overlap_len;
    //   }
    //   idxR =
    //       std::min(idxR + Window::window_overlap_len, raw_read.seq.size() -
    //       1);
    //   windows[widx].idxL = idxL;
    //   windows[widx].idxR = idxR;
    //   windows[widx].backbone_seq =
    //       std::string_view(raw_read.seq).substr(idxL, idxR - idxL + 1);
    //   if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
    //     windows[widx].backbone_qual =
    //         std::string_view(raw_read.qual).substr(idxL, idxR - idxL + 1);
    //   }
    // }
    return windows;
  }

  /**
   * @brief initialize read
   * @return void
   */
  auto make_reads() {
    reads.resize(this->raw_reads_size);
    for (auto i : std::views::iota(0ul, this->raw_reads_size)) {
      reads[i].name = std::move(raw_reads[i].name);
      reads[i].seq = std::move(raw_reads[i].seq);
      reads[i].rc_seq = std::move(raw_reads_rc_seq[i]);
      if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
        reads[i].qual = std::move(raw_reads[i].qual);
        reads[i].rc_qual = std::move(raw_reads_rev_qual[i]);
      }
      reads[i].id = name2id[reads[i].name];
      /* initial read */
    }
    // TODO: shuffle read
    // std::random_shuffle(reads.begin(), reads.end());
  }

  auto correct() {
    auto threadpool = bio::make_threadpool(threads);

    auto get_global_alignment_engine = [&](const std::size_t max_length) {
      auto aln_engine = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kNW,  // Needleman-Wunsch(global alignment)
          5,                         // match (default parameter form SPOA)
          -4,                        // mismatch
          -8,                        // gap
          -6                         // gap extension
      );
      aln_engine->Prealloc(max_length * 1.5, 5);
      return aln_engine;
    };

    auto get_local_alignment_engine = [&](const std::size_t max_length) {
      auto aln_engine = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kSW,  // Smith-Waterman(local alignment)
          5,                         // match (default parameter form SPOA)
          -4,                        // mismatch
          -8,                        // gap
          -6                         // gap extension
      );
      aln_engine->Prealloc(max_length * 1.5, 5);
      return aln_engine;
    };

    // auto init_variation_graph =
    //     [&](Window& w, std::shared_ptr<spoa::AlignmentEngine> aln_engine) {
    //   assert(aln_engine != nullptr);
    //   auto align_and_push = [&](const std::string_view seq,
    //                             const std::optional<std::string_view>& qual)
    //                             {
    //     auto aln = aln_engine->Align(seq.data(), seq.size(), w.graph);
    //     if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
    //       w.graph.AddAlignment(aln, seq.data(), seq.size(),
    //       qual.value().data(),
    //                            qual.value().size());
    //     } else {
    //       w.graph.AddAlignment(aln, seq.data(), seq.size(), 1);
    //     }
    //   };

    //   // TODO: random_view
    //   auto order = std::vector<std::size_t>(w.overlaps.size());
    //   std::iota(order.begin(), order.end(), 0);
    //   std::random_shuffle(order.begin(), order.end());
    //   align_and_push(w.backbone_seq, w.backbone_qual);
    //   for (const auto& idx : order) {
    //     align_and_push(w.overlaps[idx].t_seq, w.overlaps[idx].qual);
    //   }
    //   // ? The begin and end may need adjust due to sequencing error
    //   // ? happen at the begin and end of the read.
    //   // ? possible solution
    //   // ? 1. calcuate the coverage for first 10 and last 10 base, select
    // };

    // auto align_on_overlap_seqs =
    //     [&](Window& w, std::size_t read_idx,
    //         const std::shared_ptr<spoa::AlignmentEngine> aln_engine) {
    //   assert(aln_engine != nullptr);

    //   auto align_and_push_seq = [&](std::size_t idx, std::string_view seq,
    //                                 bool strain) {
    //     auto alignment = aln_engine->Align(seq.data(), seq.size(), w.graph);
    //     if (alignment.empty()) {
    //       return;
    //     }
    //     auto corrected = w.graph.DecodeAlignment(alignment);
    //     if (!strain) {
    //       corrected = bio::Codec::rev_comp(corrected);
    //     }
    //     std::scoped_lock lock(mutexes[idx]);
    //     // TODO: push the sequence into graph directly
    //     reads[idx].corrected_fragments.emplace_back(std::move(corrected));
    //   };
    //   for (const auto& overlap : w.overlaps) {
    //     align_and_push_seq(overlap.tid, overlap.seq, overlap.strain);
    //   }
    //   align_and_push_seq(read_idx, w.backbone_seq, true);
    // };

    auto local_aln_engines =
        std::vector<std::shared_ptr<spoa::AlignmentEngine>>{};
    auto global_aln_engines =
        std::vector<std::shared_ptr<spoa::AlignmentEngine>>{};
    for (auto i = 0ul; i < threads; i++) {
      local_aln_engines.emplace_back(
          get_local_alignment_engine(Window::window_len));
      global_aln_engines.emplace_back(
          get_global_alignment_engine(Window::window_len));
    }
    spdlog::debug("After create alignment engine");

    auto batch_job = [&](Read& read) {
      auto thread_id = threadpool.get_worker_id();
      auto local_aln_engine = local_aln_engines[thread_id];
      auto global_aln_engine = global_aln_engines[thread_id];

      // spdlog::debug("run Read {}", read.id);
      auto build_window = [&](Window& w) {
        w.build_variation_graph(global_aln_engine);
      };
      for (auto& w : read.windows) {
        build_window(w);
      }
      {
        static std::atomic_int cnt = 0;
        spdlog::debug("cnt = {}, read {} done, size = {}", ++cnt, read.id,
                      read.windows.size());
      }
    };

    // auto batch_job = [&](Read& read, std::size_t start,
    //                      std::size_t end) {
    //   auto build_window = [&](Window& w) {
    //     init_msa(w);
    //     w.graph = w.graph.PruneGraph(w.get_prune_len());
    //     align_on_overlap_seqs(w, read.idx);
    //   };
    //   for (auto i = start; i < end; i++) {
    //     build_window(read.windows[i]);
    //   }
    //   if (end == read.windows.size()) {
    //     static std::atomic_int cnt = 0;
    //     spdlog::info("cnt = {}, read {} done, size = {}", ++cnt, read.idx,
    //                 read.windows.size());
    //   }
    // };

    // auto batch_job = [&](std::size_t start, std::size_t end) {
    //   // TODO: init alignment engine here

    //   auto local_alignment_engine =
    //       get_local_alignment_engine(0);
    //   auto global_alignment_engine =
    //       get_global_alignment_engine(0);
    //   spdlog::debug("global: {}", fmt::ptr(global_alignment_engine.get()));
    //   spdlog::debug("local: {}", fmt::ptr(local_alignment_engine.get()));

    //   for (auto i = start; i < end; i++) {
    //     auto& read = reads[i];
    //     if (read.windows.empty()) {
    //       continue;
    //     }
    //     for (auto& w : read.windows) {
    //       // auto st = std::chrono::steady_clock::now();
    //       init_msa(w, global_alignment_engine);

    //       // auto subgraph = w.graph.PruneGraph(w.get_prune_len());
    //       // w.graph.Clear();
    //       // w.graph = std::move(subgraph);
    //       align_on_overlap_seqs(w, read.idx, local_alignment_engine);
    //       // auto ed = std::chrono::steady_clock::now();
    //       // if (read.idx == 0) {
    //       //   spdlog::debug("read {} window {} done, time = {}ms", read.idx,
    //       //                 w.idxL,
    //       std::chrono::duration_cast<std::chrono::milliseconds>(ed -
    //       st).count());
    //       //                 fs::path path = fmt::format(
    //       //       "/mnt/ec/ness/yolkee/thesis/tests/tmp/prune/{}",
    //       read.name);
    //       //   w.graph.PrintDot(path / fmt::format("{}-{}.dot", w.idxL,
    //       w.idxR));
    //       // }
    //     }
    //     static std::atomic_int cnt = 0;
    //     spdlog::info("cnt = {}, read {} done, size = {}", ++cnt, read.idx,
    //                  read.windows.size());
    //   }
    // };

    spdlog::info("Making windows for each read...");
    auto total_windows = 0ul;
    auto total_seqs = 0ul;
    const int take_reads = 1;  // reads.size();
#pragma omp parallel for num_threads(threads) \
    reduction(+ : total_windows, total_seqs)
    // for (int i = 0; i < take_reads; i++) {
    //   auto& read = reads[i];
    for (auto& read : reads) {
      read.windows = make_windows(read);
      total_windows += read.windows.size();
      for (auto& w : read.windows) {
        total_seqs += w.overlap_seqs.size();
      }
    }
    spdlog::debug("building window down");
    std::exit(0);

    std::vector<std::future<void>> futures;
    spdlog::info("Start building variation graph for each reads");
    for (auto& read : reads | std::views::take(take_reads)) {
      auto [_, res] = threadpool.submit(batch_job, std::ref(read));
      futures.emplace_back(std::move(res));
    }

    for (auto& f : futures) {
      f.get();
    }

    std::exit(0);

    // for (auto prev_idx = 0, idx = 0, sum = 0; idx < (int) reads.size();
    // idx++) {
    //   sum += reads[idx].windows.size();
    //   if (idx + 1 == reads.size() ||
    //       sum + reads[idx].windows.size() > WINDOW_IN_BATCH) {
    //     // spdlog::debug("submit batch job ({}, {})", prev_idx, idx + 1);
    //     auto [_, res] = threadpool.submit(batch_job, prev_idx, idx + 1);
    //     futures.emplace_back(std::move(res));
    //     sum = 0;
    //     prev_idx = idx + 1;
    //   }
    // }

    // std::vector<std::vector<std::future<void>>> futures(reads.size());
    // for (auto i : std::views::iota(0u, reads.size())) {
    //   auto& read = reads[i];
    //   // ! fix this
    //   read.windows = make_windows(read);
    //   auto& windows = read.windows;
    //   if (windows.empty()) {
    //     continue;
    //   }
    //   for (auto i = 0u; i * WINDOW_IN_BATCH < windows.size(); i++) {
    //     spdlog::debug("Read {} submit batch job {}", read.idx, i);
    //     auto [_, res] = threadpool.submit(
    //         batch_job, std::ref(read), i * WINDOW_IN_BATCH,
    //         std::min((i + 1) * WINDOW_IN_BATCH, windows.size()));
    //     futures[read.idx].emplace_back(std::move(res));
    //   }
    // }

    // for (int i = 0; i < reads.size(); i++) {
    //   auto& read = reads[i];
    //   for (auto& f : futures[i]) {
    //     f.get();
    //   }
    // }

    //   {
    //     std::size_t fragment_size_sum = 0;
    // #pragma omp parallel for reduction(+ : fragment_size_sum)
    //     for (auto& read : reads) {
    //       std::size_t len_sum = 0;
    //       for (auto& fragment : read.corrected_fragments) {
    //         len_sum += fragment.size();
    //       }
    //       fragment_size_sum += read.corrected_fragments.size();
    //       spdlog::info("read {} corrected_fragments size = {}, average len =
    //       {}, "
    //                   "total len = {}, origin read len = {}",
    //                   read.idx, read.corrected_fragments.size(),
    //                   len_sum / read.corrected_fragments.size(), len_sum,
    //                   read.seq.size());
    //     }
    //     spdlog::info("average corrected_fragments size = {}",
    //                 fragment_size_sum / reads.size());

    //     auto path =
    //     fs::path("/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/seqs");
    //     fs::create_directories(path);
    //     for (auto& read : reads | std::views::take(10)) {
    //       auto p = path / (read.name + ".txt");
    //       spdlog::info("write to {}", p.string());
    //       std::ofstream fout(p);
    //       assert(fout.is_open() && "cannot open file");
    //       for (auto& fragment : read.corrected_fragments) {
    //         fout << fragment << '\n';
    //       }
    //     }
    //   }

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

    auto corrected_read = std::vector<bio::FastaRecord<false>>{};
    return corrected_read;
  }

  FragmentedReadCorrector(
      std::vector<R>&& raw_reads, std::vector<bio::PafRecord>&& overlaps,
      const std::string& platform,
      const int thread_num = std::thread::hardware_concurrency(),
      bool debug = false)
      : BaseReadCorrector<R>(std::move(raw_reads), std::move(overlaps),
                             platform, thread_num, debug),
        mutexes(this->raw_reads_size) {
    // ! why need this pointer
    make_reads();
    assemblers.reserve(this->raw_reads_size);
    for (auto i : std::views::iota(0ul, this->raw_reads_size)) {
      assemblers.emplace_back(ReadAssembler{std::log2(reads[i].seq.size())});
    }
  }

 private:
  const auto overlap_extend_len = 50ul;

  std::vector<ReadAssembler> assemblers;
  std::vector<Read> reads;
  std::vector<std::mutex> mutexes;
};