#pragma once

#include <edlib.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <execution>
#include <future>
#include <mutex>
#include <ranges>
#include <spoa/spoa.hpp>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

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
    // const std::size_t min_overlap_length = 500ul;

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
   * @brief
   *
   *
   * @tparam R constrained by std::ranges::range
   * @param raw_read
   * @param overlap_range
   * @return std::vector<Window>
   */
  auto make_windows_for_one_read(Read& raw_read) {
    raw_read.init_windows(param.max_window_len, param.window_extend_len);
    for (const auto& overlap : raw_read.overlap_range) {
      raw_read.add_overlap_into_window(overlap, reads[overlap.t_id]);
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

#pragma omp parallel for num_threads(threads) \
    reduction(+ : total_windows, total_seqs)
    for (int i = 0; i < take_reads; i++) {
      auto& read = reads[i];
      // for (auto& read : reads) {
      make_windows_for_one_read(read);
      total_windows += read.windows.size();
      for (auto& w : read.windows) {
        total_seqs += w.overlap_seqs.size();
        for (auto& seq : w.overlap_seqs) {
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
    assert(std::ranges::is_sorted(overlaps, {}, &Overlap::q_id) &&
           "overlaps should be sorted by query id");

    // auto target = "11_8325";
    // spdlog::debug("overlap_range.size() = {}",
    // get_overlap_range(target).size());

    // auto has_pushed_read_cnt = std::atomic_int{0};
    auto cv = std::condition_variable{};

    std::function<void(Read&)> assemble_and_write = [&](Read& raw_read) {
      // spdlog::debug("Assemble read {}, len = {}", raw_read.name,
      // raw_read.len());
      // TODO: if coverage is enough, then we can assemble the read
      // remember that if we don't build the graph for this read, then there may
      // have other read that doesn't have enough data for building graph, so we
      // may need other mechanism to build the graph for this read.
      auto ready = true;
      if (!raw_read.is_finished()) {
        ready = false;
      }
      for (auto t_id : raw_read.overlap_reads_id) {
        if (!reads[t_id].is_finished()) {
          ready = false;
        }
      }
      if (ready) {
        spdlog::debug(
            "Read {} start assemble, overlap.size() = {}, window_size = {}",
            raw_read.name, raw_read.overlap_reads_id.size(),
            raw_read.windows.size());

        auto corrected_read = raw_read.get_corrected_read();
        corrected_read.name = raw_read.name;
        raw_read.clear();
        std::scoped_lock lock(corrected_reads_mutex);
        corrected_reads.emplace_back(std::move(corrected_read));
      } else {
        threadpool.submit(assemble_and_write, std::ref(raw_read));
      }
      if (corrected_reads.size() == filtered_raw_reads_size) {
        cv.notify_one();
      }
    };

    auto window_pipeline = [&](Window& window) {
      auto tid = threadpool.get_worker_id();
      thread_local auto engine =
          get_global_alignment_engine(param.max_window_len);
      window.build_variation_graph(engine);
      auto corrected_fragments = window.get_corrected_fragments(engine);
      for (auto& seq : corrected_fragments) {
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
    };

    auto read_pipeline = [&](Read& read) {
      make_windows_for_one_read(read);
      auto& windows = read.windows;
      for (auto& w : windows) {
        auto [_, res] = threadpool.submit(window_pipeline, std::ref(w));
      }
      threadpool.submit(assemble_and_write, std::ref(read));
    };

    spdlog::debug("Pushed all read into pipeline...");
    auto order = std::vector<std::size_t>(reads.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(std::execution::par, order.begin(), order.end(), [&](int a, int b) {
      return reads[a].overlap_range.size() > reads[b].overlap_range.size();
    });

    for (auto& x : order) {
      auto& read = reads[x];
      if (read.need_corrected) {
        auto [_, res] = threadpool.submit(read_pipeline, std::ref(read));
      }
    }
    spdlog::debug("Pushed read pipeline done");

    {
      spdlog::debug("Waiting all reads been assembled");
      auto mutex = std::mutex{};
      std::unique_lock lock(mutex);
      cv.wait(lock, [&]() {
        return corrected_reads.size() == filtered_raw_reads_size;
      });
      spdlog::debug("Done");
    }

    std::exit(0);

    // std::vector<std::shared_ptr<spoa::AlignmentEngine>>
    //     global_alignment_engines(threads);
    // for (auto i = 0u; i < threads; i++) {
    //   global_alignment_engines[i] =
    //       get_global_alignment_engine(param.max_window_len);
    // }

    // auto batch_job = [&](Window &window) {
    //   auto tid = threadpool.get_worker_id();
    //   auto engine = global_alignment_engines[tid];
    //   // TODO: use thread_local
    //   window.build_variation_graph(engine);
    //   auto corrected_fragments = window.get_corrected_fragments(engine);
    //   for (auto &seq : corrected_fragments) {
    //     auto lock = std::scoped_lock(mutexes[seq.read_id]);
    //     reads[seq.read_id].add_corrected_fragment(std::move(seq));
    //   }
    //   window.clear();
    //   {
    //     static std::atomic_int cnt = 0;
    //     if (cnt % 10000 == 0) {
    //       spdlog::debug("cnt = {}", cnt);
    //     }
    //     cnt++;
    //   }
    // };

    // make_windows_for_all_read();

    // std::vector<std::future<void>> futures;
    // spdlog::info("Start building variation graph for each reads");

    // auto target = "01_10002";
    // auto target_id = name2id[target];
    // spdlog::debug("target id = {}", target_id);

    // for (auto &read : reads | std::views::take(1000)) {
    //   for (auto &w : read.windows) {
    //     auto [_, res] = threadpool.submit(batch_job, std::ref(w));
    //     futures.emplace_back(std::move(res));
    //   }
    // }
    // for (auto &f : futures) {
    //   f.get();
    // }

    // auto dir = fs::path("{}/fragments", TMP_PATH);
    // fs::remove_all(dir);
    // fs::create_directories(dir);
    // for (int i = 0; i < reads.size(); i++) {
    //   spdlog::debug("{}, fragments size = {}", i,
    //   reads[i].corrected_fragments.size()); auto p = fmt::format("{}/{}.txt",
    //   dir.string(), reads[i].name); std::ofstream fout(p); for (const auto& s
    //   :reads[i].corrected_fragments) {
    //     fout << id2name[s.read_id] << '\t';
    //     fout << s.left_bound << '\t';
    //     fout << s.right_bound << '\t';
    //     fout << s.seq << '\n';
    //   }
    // }

    spdlog::debug("Start assemble reads");
    std::vector<bio::FastaRecord<false>> corrected_reads;
    auto fout = std::ofstream(TMP_PATH "/wrong2.txt");

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

  /**
   * @brief Performs read preprocessing.
   *
   * This function processes a vector of raw reads by performing various checks
   * and modifications on each read.
   *
   * @param raw_reads The vector of raw reads to be preprocessed.
   */
  auto read_preprocess(std::vector<R>& raw_reads) {
    spdlog::info("Read preprocessing...");

    auto get_overlap_range = [&](const std::string& name) {
      auto id = name2id[name];
      auto st = std::ranges::lower_bound(overlaps, id, {}, &Overlap::q_id);
      auto ed = std::ranges::upper_bound(overlaps, id, {}, &Overlap::q_id);
      return std::ranges::subrange(st, ed);
    };

    auto valid_read = [&](Read& r) {
      if (r.seq.size() < param.min_read_length) {
        return false;
      }
      if (r.overlap_range.size() == 0) {
        return false;
      }
      return true;
    };

    reads.resize(raw_reads.size());

#pragma omp parallel for
    for (auto& r : raw_reads) {
      auto read = Read(std::move(r));
      read.id = name2id[read.name];
      read.overlap_range = get_overlap_range(read.name);
      read.need_corrected = valid_read(read);
      for (const auto& overlap : read.overlap_range) {
        read.overlap_reads_id.emplace_back(overlap.t_id);
        if (!overlap.forward_strain) {
          read.create_rc();
          assert(read.create_rc());
        }
      }
      reads[read.id] = std::move(read);
    }

    filtered_raw_reads_size =
        std::ranges::count_if(reads, &Read::need_corrected);
    assert(reads.size() == raw_reads.size());
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
    auto valid_overlap = [&](const bio::PafRecord& paf) {
      auto valid_read = [&](const std::string_view name) {
        return name2id.contains(name);
      };
      if (!valid_read(paf.q_name) || !valid_read(paf.t_name)) {
        return false;
      }
      if (paf.q_name == paf.t_name) {
        return false;
      }
      // if (paf.aln_len < param.min_overlap_length) {
      //   return false;
      // }
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
            (double) std::max(q_overlap_len, t_overlap_len);
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
#pragma omp parallel for
    for (auto& paf : overlaps) {
      /* there's no std::ranges::move_if(), so sad */
      if (valid_overlap(paf)) {
        /* add reverse complement sequence and reverse quality if needed */
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

    filtered_overlaps_size = filtered_overlaps.size();
    std::sort(std::execution::par, filtered_overlaps.begin(),
              filtered_overlaps.end());
    std::swap(this->overlaps, filtered_overlaps);
  }

  auto preprocess(std::vector<R>& raw_reads,
                  std::vector<bio::PafRecord>& overlaps) {
    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();

    /* transform read name to id */
    for (auto& read : raw_reads) {
      auto id = name2id.size();
      name2id[read.name] = id;
      id2name[id] = read.name;
    }

    overlap_preprocess(overlaps);
    read_preprocess(raw_reads);
    mutexes = std::vector<std::mutex>(reads.size());
  }

  auto print_info() {
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      spdlog::info("Read type: Fastq");
    } else {
      spdlog::info("Read type: Fasta");
    }
    spdlog::info("Total raw reads: {}", unfiltered_raw_read_size);
    spdlog::info("Total overlaps: {}", unfiltered_overlap_size);
    spdlog::info("Platform: {}", platform);
    spdlog::info("Filtered raw reads: {}", filtered_raw_reads_size);
    spdlog::info("Filtered overlaps: {}", filtered_overlaps_size);
    spdlog::info("Debug mode is {}", debug ? "on" : "off");
  }

  FragmentedReadCorrector(
      std::vector<R>&& raw_reads, std::vector<bio::PafRecord>&& overlaps,
      const std::string& platform,
      const int thread_num = std::thread::hardware_concurrency(),
      bool debug = false)
      : platform(platform),
        threads(thread_num),
        threadpool(bio::make_threadpool(this->threads)),
        debug(debug) {
    preprocess(raw_reads, overlaps);
    print_info();
  }

 private:
  /* how may reads and paf records(overlap) in original data */
  std::size_t unfiltered_raw_read_size;
  std::size_t unfiltered_overlap_size;

  /* size of filtered read and filtered overlaps */
  std::size_t filtered_raw_reads_size;
  std::size_t filtered_overlaps_size;

  /* read wrapper, contains additional information */
  /* reads.size() == mutexes.size() == raw_reads_size */
  std::vector<Read> reads;
  std::vector<std::mutex> mutexes;

  /* overlap information between `raw_reads` */
  std::vector<Overlap> overlaps;

  /* mapping from read name to read id and vice versa */
  std::map<std::string, std::size_t, std::less<>> name2id;
  std::map<std::size_t, std::string> id2name;

  /* max threads and threadpool */
  std::size_t threads;
  bio::ThreadPool<> threadpool;

  /* sequencing platform of raw_reads */
  // TODO: should be enum class
  std::string platform;

  std::vector<bio::FastaRecord<false>> corrected_reads;
  std::mutex corrected_reads_mutex;

  /* debug flag */
  bool debug = false;
};