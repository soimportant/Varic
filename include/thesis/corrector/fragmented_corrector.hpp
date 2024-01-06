#pragma once

#include <execution>
#include <future>
#include <mutex>
#include <ranges>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <edlib.h>
#include <spoa/spoa.hpp>

#include "detail/window.hpp"
#include "thesis/algo/assemble/read_assembler.hpp"
#include "thesis/corrector/detail/overlap.hpp"
#include "thesis/corrector/detail/read.hpp"
#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/corrector/detail/window.hpp"
#include "thesis/format/paf.hpp"
#include "thesis/utility/threadpool/threadpool.hpp"

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

  struct Param {
    /* minimum length of read that will be corrected */
    const std::size_t min_read_length = 500ul;

    /* minimum overlap length between query read and target read */
    const std::size_t min_overlap_length = 500ul;

    /* extend length on start and end of overlap */
    const std::size_t overlap_extend_len = 25ul;

    /* longest length of a window */
    const std::size_t max_window_len = 450ul;

    /**
     * overlap length of adjanency window. Therefore, typical length of a
     * window would be (window_len + 2 * window_extend_len)
     */
    const std::size_t window_extend_len = 25ul;

  } param;

  /**
   * @brief Get the range of bio::PafRecord with the same read name.
   *
   * @param read_name name of the read
   * @return std::ranges::subrange
   */
  auto get_overlap_range(const std::string &read_name) {
    auto read_id = name2id[read_name];
    auto st = std::ranges::lower_bound(overlaps, read_id, {}, &Overlap::q_id);
    auto ed = std::ranges::upper_bound(overlaps, read_id, {}, &Overlap::q_id);
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
  auto make_windows_for_one_read(Read &raw_read) {
    
    raw_read.init_windows(param.max_window_len, param.window_extend_len);

    for (const auto &overlap : get_overlap_range(raw_read.name)) {
      raw_read.add_overlap_into_window(overlap, reads[overlap.t_id]);
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
    // return windows;
  } 

  auto ready_for_assemble(Read& raw_read) {
    // TODO: if coverage is enough, then we can assemble the read
    // remember that if we don't build the graph for this read, then there may
    // have other read that doesn't have enough data for building graph, so we
    // may need other mechanism to build the graph for this read.

    if (!raw_read.is_finished()) {
      return false;
    }
    for (auto t_id : raw_read.overlap_reads_id) {
      if (!reads[t_id].is_finished()) {
        return false;
      }
    }
    return true;
  }

  auto window_pipeline(Window& window) {
    
    auto tid = threadpool.get_worker_id();
    auto global_aln_engine = get_global_alignment_engine(param.max_window_len);
    window.build_variation_graph(global_aln_engine);
    auto corrected_fragments =
        window.get_corrected_fragments(global_aln_engine);
    for (auto &seq : corrected_fragments) {
      auto lock = std::scoped_lock(mutexes[seq.read_id]);
      reads[seq.read_id].add_corrected_fragment(std::move(seq));
    }
    window.clear();
    reads[window.read_id].finished_windows_cnt += 1;
    {
      static std::atomic_int cnt = 0;
      if (cnt % 10000 == 0) {
        spdlog::debug("cnt = {}", cnt);
      }
      cnt++;
    }
  }

  auto read_pipeline(Read& raw_read) {
    raw_read.windows = make_windows_for_one_read(raw_read);
    auto& windows = raw_read.windows;

    for (auto& w : windows) {
      auto [_, res] = threadpool.submit(window_pipeline, std::ref(w));
      futures.emplace_back(std::move(res));
    }
    
  }

  auto make_windows_for_all_read() {
    auto total_windows = 0ul;
    auto total_seqs = 0ul;
    auto total_seqs_lens = 0ul;
    // const int take_reads = reads.size();
    const int take_reads = 20;

    spdlog::debug("Taking {} reads now", take_reads);
    spdlog::info("Making windows for each reads...");

#pragma omp parallel for num_threads(threads)                                  \
    reduction(+ : total_windows, total_seqs)
    for (int i = 0; i < take_reads; i++) {
      auto &read = reads[i];
      // for (auto& read : reads) {
      make_windows_for_one_read(read);
      total_windows += read.windows.size();
      for (auto &w : read.windows) {
        total_seqs += w.overlap_seqs.size();
        for (auto &seq : w.overlap_seqs) {
          total_seqs_lens += seq.seq.size();
        }
      }
    }
    spdlog::info("building window down");
    spdlog::debug("average windows in a read = {:.2f}",
                  static_cast<double>(total_windows) / take_reads);
    spdlog::debug("average sequence inside a window = {:.2f}",
                  static_cast<double>(total_seqs) / total_windows);
    spdlog::debug("average sequence length inside a window = {:.2f}",
                  static_cast<double>(total_seqs_lens) / total_seqs);
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
          get_global_alignment_engine(param.max_window_len);
    }

    auto batch_job = [&](Window &window) {
      auto tid = threadpool.get_worker_id();
      auto engine = global_alignment_engines[tid];
      // TODO: use thread_local
      window.build_variation_graph(engine);
      auto corrected_fragments = window.get_corrected_fragments(engine);
      for (auto &seq : corrected_fragments) {
        auto lock = std::scoped_lock(mutexes[seq.read_id]);
        reads[seq.read_id].add_corrected_fragment(std::move(seq));
      }
      window.clear();
      {
        static std::atomic_int cnt = 0;
        if (cnt % 10000 == 0) {
          spdlog::debug("cnt = {}", cnt);
        }
        cnt++;
      }
    };

    make_windows_for_all_read();

    std::vector<std::future<void>> futures;
    spdlog::info("Start building variation graph for each reads");

    auto target = "01_10002";
    auto target_id = name2id[target];
    spdlog::debug("target id = {}", target_id);

    for (auto &read : reads | std::views::take(20)) {
      for (auto &w : read.windows) {
        auto [_, res] = threadpool.submit(batch_job, std::ref(w));
        futures.emplace_back(std::move(res));
      }
    }
    for (auto &f : futures) {
      f.get();
    }

    // auto dir = fs::path("{}/fragments", TMP_PATH);
    // fs::remove_all(dir);
    // fs::create_directories(dir);
    // for (int i = 0; i < reads.size(); i++) {
    //   spdlog::debug("{}, fragments size = {}", i, reads[i].corrected_fragments.size()); 
    //   auto p = fmt::format("{}/{}.txt", dir.string(), reads[i].name);
    //   std::ofstream fout(p);
    //   for (const auto& s :reads[i].corrected_fragments) {
    //     fout << id2name[s.read_id] << '\t';
    //     fout << s.left_bound << '\t';
    //     fout << s.right_bound << '\t';
    //     fout << s.seq << '\n';
    //   }
    // }

    spdlog::debug("Start assemble reads");
    std::vector<bio::FastaRecord<false>> corrected_reads;
    auto fout = std::ofstream(TMP_PATH"/wrong2.txt");
    

    
// #pragma omp parallel for num_threads(threads)
//     for (auto &read : reads) {
//       // spdlog::debug("Assemble read {}, len = {}", id2name[read.id],
//       // read.len());
//       auto corrected_read = read.get_corrected_read();
//       corrected_read.name = id2name[read.id];

//       #pragma omp critical
//       if (corrected_read.seq.size() != 0) {
//         corrected_reads.emplace_back(std::move(corrected_read));
//         if (corrected_reads.size() % 100 == 0) {
//           spdlog::debug("Corrected {} reads", corrected_reads.size());
//         }
//       }
//     }

    return corrected_reads;

    // fs::remove_all("/mnt/ec/ness/yolkee/thesis/tests/fragments");
    // fs::create_directories("/mnt/ec/ness/yolkee/thesis/tests/fragments");
    // for (int i = 0; i < reads.size(); i++) {
    //   spdlog::debug("{}, fragments size = {}", i,
    //   reads[i].corrected_fragments.size()); auto p =
    //   fmt::format("/mnt/ec/ness/yolkee/thesis/tests/fragments/{}.txt",reads[i].name);
    //   std::ofstream fout(p);
    //   for (const auto& s :reads[i].corrected_fragments) {
    //     fout << s.left_bound << '\t' << s.right_bound << '\t' << s.seq <<
    //     '\n';
    //   }
    // }

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
      std::vector<R> &&raw_reads, std::vector<bio::PafRecord> &&overlaps,
      const std::string &platform,
      const int thread_num = std::thread::hardware_concurrency(),
      bool debug = false)
      : platform(platform), threads(thread_num),
        threadpool(bio::make_threadpool(this->threads)), debug(debug) {
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
  auto read_preprocess(std::vector<R> &raw_reads) {
    spdlog::info("Read preprocessing...");

    /* there's no std::ranges::move_if() */
    for (auto &r : raw_reads) {
      if (r.seq.size() >= param.min_read_length) {
        reads.emplace_back(std::move(r));
      }
    }

    // compile error
    // std::sort(std::execution::par, reads.begin(), reads.end());
    raw_reads_size = reads.size();

    /* transform read name to id */
    for (auto &raw_read : reads) {
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
  auto overlap_preprocess(std::vector<bio::PafRecord> &overlaps) {
    spdlog::info("Overlap preprocessing...");

    // In original .paf file of overlaps between raw_reads, the match base
    // devided by alignment length is very low, figure out how it be computed
    // and why it is so low.
    // const auto min_overlap_identity = 0.7l;

    auto valid_overlap = [&](const bio::PafRecord &paf) {
      auto valid_read = [&](const std::string_view name) {
        return name2id.contains(name);
      };
      if (!valid_read(paf.q_name) || !valid_read(paf.t_name)) {
        return false;
      }
      if (paf.q_name == paf.t_name) {
        return false;
      }
      if (paf.aln_len < param.min_overlap_length) {
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

        auto q_overlap_len = paf.q_end - paf.q_start + 1;
        auto t_overlap_len = paf.t_end - paf.t_start + 1;
        auto possible_match_ratio =
            std::min(q_overlap_len, t_overlap_len) /
            (double)std::max(q_overlap_len, t_overlap_len);
        if (possible_match_ratio < 0.7) {
          return false;
        }

        const double internal_threshold = 0.8;
        auto q_head_len = paf.q_start;
        auto q_tail_len = paf.q_len - paf.q_end;
        auto t_head_len = paf.t_start;
        auto t_tail_len = paf.t_len - paf.t_end;
        auto overhang_len =
            paf.strand == '+'
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

    auto filtered_overlaps = std::vector<Overlap>{};
    auto reads_need_create_rc = std::vector<int>(raw_reads_size, false);
#pragma omp parallel for
    for (auto &paf : overlaps) {
      /* there's no std::ranges::move_if(), so sad */
      if (valid_overlap(paf)) {
        /* add reverse complement sequence and reverse quality if needed */
        if (paf.strand == '-') {
          auto tid = name2id[paf.t_name];
          reads_need_create_rc[tid] = true;
        }
        auto overlap = Overlap{.q_id = name2id[paf.q_name],
                               .q_idx_L = paf.q_start,
                               .q_idx_R = paf.q_end,
                               .q_seq_len = paf.q_len,
                               .t_id = name2id[paf.t_name],
                               .t_idx_L = paf.t_start,
                               .t_idx_R = paf.t_end,
                               .t_seq_len = paf.t_len,
                               .forward_strain = (paf.strand == '+')};
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


  /* how may reads in original data */
  std::size_t unfiltered_raw_read_size;
  /* how many paf records(overlap) in original data */
  std::size_t unfiltered_overlap_size;

  /* overlap information between `raw_reads` */
  std::vector<Overlap> overlaps;

  /* size of filtered read */
  std::size_t raw_reads_size;

  /* size of filtered overlaps */
  std::size_t overlaps_size;

  /* sequencing platform of raw_reads */
  // TODO: should be enum class
  std::string platform;

  /* mapping from read name to read id and vice versa */
  std::map<std::string, std::size_t, std::less<>> name2id;
  std::map<std::size_t, std::string> id2name;

  /* threads */
  std::size_t threads;
  bio::ThreadPool<> threadpool;
  std::vector<std::future<void>> futures;


  /* read wrapper, contains additional information */
  /* reads.size() == mutexes.size() == raw_reads_size */
  std::vector<Read> reads;
  std::vector<std::mutex> mutexes;

  /* debug flag */
  bool debug = false;
};