#pragma once


#include <concepts>
#include <execution>
#include <future>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <ranges>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <spoa/spoa.hpp>
#include <edlib.h>

#include "thesis/algo/assemble/read_assembler.hpp"
#include "thesis/corrector/detail/read.hpp"
#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/corrector/detail/window.hpp"
#include "thesis/format/paf.hpp"
#include "thesis/utility/threadpool/threadpool.hpp"

/**
 * @brief align overlap part between query read and target read, and return the
 * cigar string
 *
 * @param q subsequence on query read
 * @param t subsequence on target read
 * @return bio::Cigar
 */
auto align(const std::string_view q, const std::string_view t) {
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
                                 const bio::PafRecord& overlap) {
  auto window_cnt = std::ranges::size(window_range);
  auto [q_start, q_end] = std::make_pair(overlap.q_start, overlap.q_end);
  auto [t_start, t_end] = std::make_pair(overlap.t_start, overlap.t_end);

  /* corresponding interval of target read and query read in each window */
  auto t_breakpoints =
      std::vector<std::pair<std::size_t, std::size_t>>(window_cnt);
  auto q_breakpoints =
      std::vector<std::pair<std::size_t, std::size_t>>{window_cnt};

  /* iterator for query read and target read */
  auto [q_iter, t_iter] = std::make_pair(q_start, t_start);
  auto q_last_match = q_iter;
  auto t_last_match = t_iter;

  /* store the window indexes that window is not begin at 'M' region */
  auto wait_for_first_match = std::queue<std::size_t>{};

  /* store the window indexes that not reach the end of window */
  auto wait_for_last_match = std::queue<std::size_t>{};

  /**
   * @brief set the start position of target read when meet a match or mismatch
   * in cigar string for each window.
   */
  auto set_start_pos_when_match = [&](std::size_t idx) {
    t_breakpoints[idx].first = t_iter;
    q_breakpoints[idx].first = q_iter;
    wait_for_last_match.push(idx);
  };

  /**
   * @brief when q_iter increase, check whether the window inside
   * `wait_for_first_match` and `wait_for_last_match` reach the end or not.
   *
   * for `wait_for_first_match`, if window not begin at 'M' region but reach
   * the end, set interval as an empty sequence.
   *
   * for `wait_for_last_match`, if window reach the end, record the last match
   * position of query read and target read.
   */
  auto check_whether_reach_end = [&]() {
    while (!wait_for_first_match.empty()) {
      auto idx = wait_for_first_match.front();
      if (q_iter + 1 == std::min(q_end, window_range[idx].end)) {
        wait_for_first_match.pop();
        q_breakpoints[idx].first = q_breakpoints[idx].second = q_iter + 1;
      } else {
        break;
      }
    }
    while (!wait_for_last_match.empty()) {
      auto idx = wait_for_last_match.front();
      /* q_end < window_range[idx].end -> the last window */
      if (q_iter + 1 == std::min(q_end, window_range[idx].end)) {
        wait_for_last_match.pop();
        t_breakpoints[idx].second = t_last_match + 1;
        q_breakpoints[idx].second = q_last_match + 1;
      } else {
        break;
      }
    }
  };

  auto w_idx = 0u;
  for (const auto& [len, op] : cigar) {
    if (op == 'M') {
      while (!wait_for_first_match.empty()) {
        auto idx = wait_for_first_match.front();
        wait_for_first_match.pop();
        set_start_pos_when_match(idx);
      }
      for (auto i = 0; i < len; i++, q_iter++, t_iter++) {
        /**
         * q_start > window_range[w_idx].start -> the first window
         *
         * while loop may run twice when two window start at same position,
         * which means the length of first window is less than
         * `Window::window_extend_len`
         */
        q_last_match = q_iter;
        t_last_match = t_iter;
        while (w_idx < window_cnt &&
               q_iter == std::max(q_start, window_range[w_idx].start)) {
          set_start_pos_when_match(w_idx);
          w_idx++;
        }
        check_whether_reach_end();
      }
    } else if (op == 'I') {
      for (auto i = 0; i < len; i++, q_iter++) {
        while (w_idx < window_cnt &&
               q_iter == std::max(q_start, window_range[w_idx].start)) {
          wait_for_first_match.push(w_idx);
          w_idx++;
        }
        check_whether_reach_end();
      }
    } else if (op == 'D') {
      t_iter += len;
    } else {
      spdlog::error("unknown cigar operation: {}", op);
      std::exit(EXIT_FAILURE);
    }
  }

  /* transform coordinate system to target read if relative is reverse */
  if (overlap.strand == '-') {
    for (auto i = 0ul; i < window_cnt; i++) {
      auto len = t_breakpoints[i].second - t_breakpoints[i].first;
      auto offset = t_breakpoints[i].first - t_start;
      t_breakpoints[i] = std::make_pair(t_end - offset - len, t_end - offset);
    }
  }

  /* transform coordinate system relative to each window */
  for (auto i = 0; i < window_cnt; i++) {
    // spdlog::debug("({}, {}) <-> ({}, {})", q_breakpoints[i].first,
    //               q_breakpoints[i].second, window_range[i].start,
    //               window_range[i].end);
    assert(q_breakpoints[i].first >= window_range[i].start &&
           q_breakpoints[i].first <= window_range[i].end &&
           "q_breakpoints[i].first out of range");
    q_breakpoints[i].first -= window_range[i].start;
    q_breakpoints[i].second -= window_range[i].start;
  }

  return std::make_pair(std::move(q_breakpoints), std::move(t_breakpoints));
}

template <class R>
  requires std::derived_from<R, bio::FastaRecord<R::encoded>>
class FragmentedReadCorrector {
 private:
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

    /* initialize boundaries and backbone sequence for each window */
    for (auto pos = 0ul; pos < raw_read.seq.size(); pos += param.window_len) {
      auto w_idx = windows.size();
      auto window = Window{raw_read.id, w_idx, pos, pos + param.window_len};
      window.extend(param.window_extend_len, raw_read.seq.size());

      window.backbone.read_id = raw_read.id;
      window.backbone.left_bound = window.start;
      window.backbone.right_bound = window.end;
      std::tie(window.backbone.seq, window.backbone.qual) =
          raw_read.subview(true, window.start, window.end);
      window.backbone.forward_strain = true;

      // ? if the last window is too short, that means the right bound of the
      // ? second to last window is the end of read, then should we remove the
      // ? last window?
      windows.emplace_back(std::move(window));
    }

    auto get_window_range = [&](std::size_t start, std::size_t end) {
      auto st = std::ranges::upper_bound(windows, start, {}, &Window::start);
      if (st != windows.end()) {
        st = std::prev(st);
      }
      auto ed = std::ranges::lower_bound(windows, end, {}, &Window::start);
      return std::ranges::subrange(st, ed);
    };

    for (const auto& overlap : get_overlap_range(raw_read.name)) {
      auto t_id = name2id[overlap.t_name];
      auto forward_strain = (overlap.strand == '+');

      auto [q_start, q_end] = std::make_pair(overlap.q_start, overlap.q_end);
      auto [t_start, t_end] = std::make_pair(overlap.t_start, overlap.t_end);
      auto [q_seq_view, q_qual_view] = raw_read.subview(true, q_start, q_end);
      auto [t_seq_view, t_qual_view] =
          reads[t_id].subview(forward_strain, t_start, t_end);
      auto window_range = get_window_range(q_start, q_end);
      auto window_idx = window_range.begin() - windows.begin();

      auto cigar = align(q_seq_view, t_seq_view);
      auto [q_breakpoints, t_breakpoints] =
          find_breakpoints_from_cigar(cigar, window_range, overlap);
      const auto window_cnt = std::ranges::size(window_range);
      assert(q_breakpoints.size() == t_breakpoints.size() &&
             "size not equal between q_breakpoints and t_breakpoints");

      for (auto i = 0u; i < window_cnt; i++, window_idx++) {
        auto& window = windows[window_idx];
        auto& [t_window_st, t_window_ed] = t_breakpoints[i];
        if (t_window_st >= t_end) {
          continue;
        }
        assert(t_window_ed >= t_window_st &&
               "t_window_ed < t_window_st, wrong boundary");
        auto [t_window_seq, t_window_qual] =
            reads[t_id].subview(forward_strain, t_window_st, t_window_ed);

        // when racon add sequence into window, it will store the pos of first
        // match, then sort the overlap_seqs according the pos of first match
        window.add_sequence(Sequence<>{.read_id = t_id,
                                       .left_bound = t_window_st,
                                       .right_bound = t_window_ed,
                                       .seq = t_window_seq,
                                       .qual = t_window_qual,
                                       .forward_strain = forward_strain},
                            q_breakpoints[i]);
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


 public:
  /**
   * @brief correct the read
   *
   * @return auto
   */
  auto correct() {
    auto threadpool = bio::make_threadpool(threads);



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

    std::vector<std::shared_ptr<spoa::AlignmentEngine>>
        global_alignment_engines(threads);
    for (auto i = 0u; i < threads; i++) {
      global_alignment_engines[i] =
          get_global_alignment_engine(param.window_len);
    }

    // auto batch_job = [&](Read& read) {
    //   auto tid = threadpool.get_worker_id();
    //   // auto local_aln_engine = local_alignment_engines[tid];
    //   auto global_aln_engine = global_alignment_engines[tid];

    //   auto build_window = [&](Window& w) {
    //     w.build_variation_graph(global_aln_engine);

    //     // auto corrected_fragments = w.get_corrected_fragments(global_aln_engine);
    //     // for (auto& sequence : w.overlap_seqs) {
    //     //   auto t_id = sequence.read_id;
    //     //   // TODO: check the siz e of all pruned graph size
    //     //   // TODO: is it small than window_len?
    //     //   // TODO: or there's a lot of graph is not be pruned?
    //     //   auto fragment =
    //     //       w.get_corrected_fragment(local_aln_engine, sequence.seq);
    //     //   // std::scoped_lock lock(mutexes[t_id]);
    //     //   // reads[t_id].corrected_fragments.emplace_back(std::move(fragment));
    //     // }
    //   };
    //   for (auto& w : read.windows) {
    //     build_window(w);
    //   }
     
    // };

    auto batch_job = [&](Window& window) {
      auto tid = threadpool.get_worker_id();
      // auto local_aln_engine = local_alignment_engines[tid];
      auto global_aln_engine = global_alignment_engines[tid];

      window.build_variation_graph(global_aln_engine);

      auto corrected_fragments = window.get_corrected_fragments(global_aln_engine);
      for (auto& seq : corrected_fragments) {
        auto lock = std::scoped_lock(mutexes[seq.read_id]);
        reads[seq.read_id].add_corrected_fragment(seq);
      }
      corrected_fragments.clear();
      window.clear();
      {
        static std::atomic_int cnt = 0;
        if (cnt % 10000 == 0) {
          spdlog::debug("cnt = {}", cnt);
        }
        cnt++;
      }
    };

    auto total_windows = 0ul;
    auto total_seqs = 0ul;
    auto total_seqs_lens = 0ul;
    const int take_reads = reads.size();

    spdlog::debug("Taking {} reads now", take_reads);
    spdlog::info("Making windows for each reads...");

#pragma omp parallel for num_threads(threads) \
    reduction(+ : total_windows, total_seqs)
    for (int i = 0; i < take_reads; i++) {
      auto& read = reads[i];
      // for (auto& read : reads) {
      read.windows = make_windows(read);
      total_windows += read.windows.size();
      for (auto& w : read.windows) {
        total_seqs += w.overlap_seqs.size();
        for (auto& seq : w.overlap_seqs) {
          total_seqs_lens += seq.seq.size();
        }
      }
    }
    spdlog::info("building window down");
    spdlog::debug("average windows in a read = {}", total_windows / take_reads);
    spdlog::debug("average sequence inside a window = {:.2f}",
                  static_cast<double>(total_seqs) / total_windows);
    spdlog::debug("average sequence length inside a window = {:.2f}",
                  static_cast<double>(total_seqs_lens) / total_seqs);

    std::vector<std::future<void>> futures;

    spdlog::info("Start building variation graph for each reads");
    for (auto& read : reads | std::views::take(take_reads)) {
      for (auto&& window : read.windows) {
        auto [_, res] = threadpool.submit(batch_job, std::ref(window));
        futures.emplace_back(std::move(res));
      }
    }
    for (auto& f : futures) {
      f.get();
    }





    
    // {
    //   auto nodes_fout = std::ofstream("/mnt/ec/ness/yolkee/thesis/tests/stat/nodes.txt");
    //   auto edges_fout = std::ofstream("/mnt/ec/ness/yolkee/thesis/tests/stat/edges.txt");
    //   for (auto i = 0u; i < reads.size(); i++) {
    //     for (auto& w : reads[i].windows) {
    //       auto [nodes_size, edges_size] = w.get_graph_info();
    //       nodes_fout << nodes_size << ' ';
    //       edges_fout << edges_size << ' ';
    //     }
    //   }
    //   nodes_fout << '\n';
    //   edges_fout << '\n';
    // }

    // {
    //   auto fout = std::ofstream("/mnt/ec/ness/yolkee/thesis/tests/stat/graph.txt");
    //   for (auto i = 0; i < reads.size(); i++) {
    //     auto& read = reads[i];
    //     fout << read.len() << ' ' << read.assembler.get_graph_size() << '\n';
    //   }
    // }

    // for (int i = 0; i < reads.size(); i++) {
    //   spdlog::debug("{}, fragments size = {}", i, reads[i].corrected_fragments.size()); 
    //   auto p = fmt::format("/mnt/ec/ness/yolkee/thesis/tests/fragments/{}.txt",reads[i].name); 
    //   std::ofstream fout(p); 
    //   for (const auto& s :reads[i].corrected_fragments) {
    //     fout << s.left_bound << '\t' << s.right_bound << '\t' << s.seq << '\n';
    //   }
    // }
    

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

  auto print_info() -> void {
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      spdlog::info("Read type: Fastq");
    } else {
      spdlog::info("Read type: Fasta");
    }
    spdlog::info("Total raw reads: {}", unfiltered_raw_read_size);
    spdlog::info("Total overlaps: {}", unfiltered_overlap_size);
    spdlog::info("Platform: {}", platform);
    spdlog::info("Filtered raw reads: {}", raw_reads_size);
    spdlog::info("Filtered overlaps: {}", overlaps_size);
    spdlog::info("Debug mode is {}", debug ? "on" : "off");
  }

  FragmentedReadCorrector(
      std::vector<R>&& raw_reads, std::vector<bio::PafRecord>&& overlaps,
      const std::string& platform,
      const int thread_num = std::thread::hardware_concurrency(),
      bool debug = false)
      : platform(platform), threads(thread_num), debug(debug) {
    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();

    read_preprocess(raw_reads);
    overlap_preprocess(overlaps);
    print_info();

    mutexes = std::vector<std::mutex>(raw_reads_size);
  }

 private:
  /**
   * @brief
   *
   * @param reads
   * @return auto
   */
  auto read_preprocess(std::vector<R>& raw_reads) {
    spdlog::info("Read preprocessing...");

    /* there's no std::ranges::move_if() */
    for (auto& r : raw_reads) {
      if (r.seq.size() >= param.min_read_length) {
        reads.emplace_back(std::move(r));
      }
    }
    // TODO: sort by name is not a proper way, we may face some bias at assembly
    // stage
    std::sort(std::execution::par, reads.begin(), reads.end());
    raw_reads_size = reads.size();

    /* transform read name to id */
    for (auto& raw_read : reads) {
      raw_read.id = name2id.size();
      name2id[raw_read.name] = raw_read.id;
      id2name[raw_read.id] = raw_read.name;
    }
  }

  /**
   * @brief
   * @param raw_reads
   * @param overlaps
   * @return
   */
  auto overlap_preprocess(std::vector<bio::PafRecord>& overlaps) {
    spdlog::info("Overlap preprocessing...");

    // In original .paf file of overlaps between raw_reads, the match base
    // devided by alignment length is very low, figure out how it be computed
    // and why it is so low.
    // const auto min_overlap_identity = 0.7l;

    auto valid_overlap = [&](const bio::PafRecord& overlap) {
      auto valid_read = [&](const std::string_view name) {
        return std::ranges::binary_search(reads, name, {}, &Read::name);
      };
      if (!valid_read(overlap.q_name) || !valid_read(overlap.t_name)) {
        return false;
      }
      if (overlap.q_name == overlap.t_name) {
        return false;
      }
      if (overlap.aln_len < param.min_overlap_length) {
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
    auto reads_need_create_rc = std::vector<int>(raw_reads_size, false);
    /* there's no std::ranges::move_if(), so sad */
#pragma omp parallel for
    for (auto& overlap : overlaps) {
      if (valid_overlap(overlap)) {
        /* add reverse complement sequence and reverse quality if needed */
        if (overlap.strand == '-') {
          auto tid = name2id[overlap.t_name];
          reads_need_create_rc[tid] = true;
        }
        overlap.extend(param.overlap_extend_len);
#pragma omp critical
        filtered_overlaps.emplace_back(std::move(overlap));
      }
    }

#pragma omp parallel for
    for (auto i : std::views::iota(0ul, raw_reads_size)) {
      if (reads_need_create_rc[i]) {
        reads[i].create_rc();
      }
    }
    overlaps_size = filtered_overlaps.size();
    std::sort(std::execution::par, filtered_overlaps.begin(),
              filtered_overlaps.end());
    std::swap(this->overlaps, filtered_overlaps);
  }

  struct Param {
    /* minimum length of read that will be corrected */
    const std::size_t min_read_length = 500ul;

    /* minimum overlap length between query read and target read */
    const std::size_t min_overlap_length = 500ul;

    /* extend length on start and end of overlap */
    const std::size_t overlap_extend_len = 50ul;

    /* length of a window */
    const std::size_t window_len = 470ul;

    /**
     * overlap length of adjanency window. Therefore, typical length of a
     * window would be (window_len + 2 * window_extend_len)
     */
    const std::size_t window_extend_len = 15ul;

  } param;

  /* how much reads in original data */
  std::size_t unfiltered_raw_read_size;
  /* how much overlap records in original data */
  std::size_t unfiltered_overlap_size;

  /* overlap information between `raw_reads` */
  std::vector<bio::PafRecord> overlaps;

  /* size of filtered read */
  std::size_t raw_reads_size;
  /* size of filtered overlaps */
  std::size_t overlaps_size;

  /* sequencing platform of raw_reads */
  // TODO: should be enum class
  std::string platform;

  /* mapping from read name to read id and vice versa */
  std::map<std::string, std::size_t> name2id;
  std::map<std::size_t, std::string> id2name;

  /* threads */
  std::size_t threads;

  /* debug flag */
  bool debug = false;

  /* read wrapper, contains additional information */
  /* reads.size() == mutexes.size() == raw_reads_size */
  std::vector<Read> reads;
  std::vector<std::mutex> mutexes;
};