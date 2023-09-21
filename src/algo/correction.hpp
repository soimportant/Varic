#pragma once

#include <algorithm>
#include <concepts>
#include <execution>
#include <map>
#include <mutex>
#include <random>
#include <ranges>
#include <string>
#include <thread>
#include <vector>

#include <biovoltron/file_io/all.hpp>
#include <spdlog/spdlog.h>
#include <spoa/spoa.hpp>

// ! ugly solution, fix it later
#include "../utility/all.hpp"
#include "assembler.hpp"

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

    const auto min_overlap_length = 1000ull;
    // In original .paf file of overlaps between raw_reads, the match base
    // devided by alignment length is very low, figure out how it be computed
    // and why it is so low.
    // const auto min_overlap_identity = 0.7l;
    auto valid_overlap = [&](const bio::PafRecord& overlap) {
      auto valid_read = [&](const std::string_view name) {
        return std::ranges::binary_search(raw_reads, name, {},
                                          &bio::FastaRecord<false>::name);
      };
      if (!valid_read(overlap.qname) || !valid_read(overlap.tname)) {
        return false;
      }
      if (overlap.qname == overlap.tname) {
        return false;
      }
      if (overlap.aln_len < min_overlap_length) {
        return false;
      }
      // if (overlap.aln_len * min_overlap_identity > overlap.match) {
      //   return false;
      // }
      return true;
    };

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
    Window() = default;

    // delete copy constructor and copy assignment
    Window(const Window&) = delete;
    Window& operator=(const Window&) = delete;

    // default move constructor and move assignment
    Window(Window&&) = default;
    Window& operator=(Window&&) = default;

    // basic length of a window
    static const auto window_len = 500ul;
    // overlap length of adjanency window. Therefore, typeical length of a
    // single window should be window_len + 2 * window_overlap_len
    static const auto window_overlap_len = 50ul;
    // extend length of overlap
    static const auto overlap_extend_len = 25ul;
    // minimum length of overlap
    static const auto overlap_min_len = window_len - window_overlap_len;

    /* backbone read information */
    std::size_t idxL, idxR;
    std::string backbone_seq;

    /* overlap information of backbone read */
    std::vector<std::size_t> overlap_read_idx;
    std::vector<std::pair<std::size_t, std::size_t>> overlap_boundaries;
    std::vector<std::string> overlap_seqs;
    std::vector<bool> overlap_strains;

    spoa::Graph graph;

    auto get_prune_len() const noexcept -> std::size_t {
      auto len = backbone_seq.size();
      return len;
    }

    auto print() const noexcept {
      spdlog::debug("idxL = {}, idxR = {}", idxL, idxR);
      spdlog::debug("backbone_seq.size() = {}", backbone_seq.size());
      auto sum = std::accumulate(overlap_seqs.begin(), overlap_seqs.end(), 0ul,
                                 [](auto a, auto& b) { return a + b.size(); });
      spdlog::debug("average overlap_seqs len = {}", sum / overlap_seqs.size());
      spdlog::debug("get prune len = {}", get_prune_len());
    }
  };

  class Read : public bio::FastqRecord<false> {
  public:
    Read() = default;

    // disable copy constructor and copy assignment
    Read(const Read&) = delete;
    Read& operator=(const Read&) = delete;

    // default move constructor and move assignment
    Read(Read&&) = default;
    Read& operator=(Read&&) = default;

    // TODO: coverage -> segment tree like data structure
    // ? we may add a data structure here for recording the coverage covered
    // ? by corrected_fragments, if the coverage is enough, we can assemble
    // ? the corrected read without building MSA.

    std::size_t idx;

    // reverse complement sequence
    std::string rc_seq;

    // correct seqs from windows of others reads, may further assemble to
    // correct version of this read
    std::vector<std::string> corrected_fragments;

    // windows of this read
    std::vector<Window> windows;
  };

  /**
   * @brief Get the overlap range object
   *
   * @param read_name
   * @return std::range::range
   */
  std::ranges::range auto get_overlap_range(const std::string_view read_name) {
    auto st = std::ranges::lower_bound(overlaps, read_name, {},
                                       &bio::PafRecord::qname);
    auto ed = std::ranges::upper_bound(overlaps, read_name, {},
                                       &bio::PafRecord::qname);
    return std::ranges::subrange(st, ed);
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
      /* boundaries information of target read -> (widx, (lbound, rbound)) */
      auto boundaries = std::vector<
          std::pair<std::size_t, std::pair<std::size_t, std::size_t>>>{};
      auto tid = name2id[overlap.tname];

      /* extend overlap boundaries at query read and target read */
      auto qstart = overlap.qstart > Window::overlap_extend_len
                        ? overlap.qstart - Window::overlap_extend_len
                        : 0ul;
      auto qend =
          std::min(overlap.qlen - 1, overlap.qend + Window::overlap_extend_len);

      auto tstart = overlap.tstart > Window::overlap_extend_len
                        ? overlap.tstart - Window::overlap_extend_len
                        : 0ul;
      auto tend =
          std::min(overlap.tlen - 1, overlap.tend + Window::overlap_extend_len);

      /* initialize boundaries */
      auto init_boundaries = [&]() {
        for (auto widx = qstart / Window::window_len;
             widx * Window::window_len < qend; ++widx) {
          auto qidxL = std::max(qstart, widx * Window::window_len);
          auto qidxR = std::min((widx + 1) * Window::window_len - 1, qend);
          auto tidxL = std::min(tstart + (qidxL - qstart), tend);
          auto tidxR = std::min(tstart + (qidxR - qstart), tend);
          if (tidxL == tidxR) {
            continue;
          }
          boundaries.emplace_back(widx, std::make_pair(tidxL, tidxR));
        }
      };

      /**
       * adjust boundaries according to length difference between query
       * sequence and target sequence
       */
      auto adjust_boundaries = [&]() {
        auto len_diff = (std::int64_t)(tend - tstart + 1) -
                        (std::int64_t)(qend - qstart + 1);
        auto offset = len_diff / (std::int64_t) boundaries.size();
        for (auto i = 1u; i < boundaries.size(); i++) {
          auto& prev_tidxR = boundaries[i - 1].second.second;
          auto& curr_tidxL = boundaries[i].second.first;

          // do this way for avoiding unsigned integer underflow
          // I don't want to debug...
          if (offset >= 0) {
            prev_tidxR = std::min(prev_tidxR + offset * i, tend);
            curr_tidxL = std::min(curr_tidxL + offset * i, tend);
          } else {
            auto n_offset = -offset;
            prev_tidxR =
                (prev_tidxR > n_offset * i) ? prev_tidxR - n_offset * i : 0ul;
            curr_tidxL =
                (curr_tidxL > n_offset * i) ? curr_tidxL - n_offset * i : 0ul;
          }
        }
        auto& [last_tidxL, last_tidxR] = boundaries.back().second;
        last_tidxR = std::min(last_tidxR + len_diff, tend);
        if (last_tidxR <= last_tidxL) {
          boundaries.pop_back();
        }
      };

      /* extend boundaries for making adjanency window have overlap */
      auto extend_boundaries = [&]() {
        std::ranges::for_each(boundaries, [&](auto& boundary) {
          auto& [_, tidx] = boundary;
          auto& [tidxL, tidxR] = tidx;
          /* cannot use std::min due to unsigned integer subtraction */
          tidxL = tidxL > Window::window_overlap_len
                      ? tidxL - Window::window_overlap_len
                      : 0ul;
          tidxR = std::min(tidxR + Window::window_overlap_len, tend);
        });
      };

      auto print = [&]() {
        std::stringstream ss;
        ss << overlap;
        spdlog::debug("{}", ss.str());
        auto len_diff = (std::int64_t)(tend - tstart + 1) -
                        (std::int64_t)(qend - qstart + 1);
        auto offset = len_diff / (std::int64_t) boundaries.size();
        spdlog::debug("len_diff = {}, offset = {}", len_diff, offset);
        for (auto& [widx, tidx] : boundaries) {
          spdlog::debug("widx = {}, tidxL = {}, tidxR = {}", widx,
          tidx.first,
                        tidx.second);
        }
      };

      /* push sequence into window */
      auto push_sequence = [&]() {
        for (const auto& [widx, tidx] : boundaries) {
          while (widx >= windows.size()) {
            windows.emplace_back(Window{});
          }
          auto forward_strain = (overlap.strand == '+');
          auto tseq = std::string_view{};
          tseq = forward_strain ? reads[tid].seq : reads[tid].rc_seq;
          auto [tidxL, tidxR] = tidx;
          if (tidxL > tidxR) {
            print();
          }
          assert(tidxL <= tidxR && "wrong boundary");
          assert(tidxL < tend && "tidxL out of range");

          auto subseq = tseq.substr(tidxL, tidxR - tidxL + 1);
          assert(subseq.size() < Window::window_len * 2 &&
                 "overlap sequence too long");
          if (subseq.size() < Window::overlap_min_len) {
            continue;
          }
          // TODO: Directly push sequence into graph
          windows[widx].overlap_read_idx.emplace_back(tid);
          windows[widx].overlap_boundaries.emplace_back(tidxL, tidxR);
          windows[widx].overlap_seqs.emplace_back(std::move(subseq));
          windows[widx].overlap_strains.emplace_back(forward_strain);
        }
      };

      init_boundaries();
      adjust_boundaries();
      extend_boundaries();
      push_sequence();
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

  /**
   * @brief initialize read
   * @return void
   */
  auto make_reads() {
    reads.resize(raw_read_size);
    for (auto i : std::views::iota(0ul, raw_read_size)) {
      reads[i].name = std::move(raw_reads[i].name);
      reads[i].seq = std::move(raw_reads[i].seq);
      reads[i].rc_seq = std::move(raw_reads_rc_seq[i]);
      reads[i].idx = name2id[reads[i].name];
      /* initial read */
    }
    // TODO: shuffle read
    // std::random_shuffle(reads.begin(), reads.end());
  }


  auto correct() {
    constexpr auto WINDOW_IN_BATCH = 1500ul;
    auto threadpool = bio::make_threadpool(threads);

    auto get_global_alignment_engine = [&](const std::size_t max_length) {
      auto engine = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kNW, // Needleman-Wunsch(global alignment)
          5,                        // match (default parameter form SPOA)
          -4,                       // mismatch
          -8,                       // gap
          -6                        // gap extension
      );
      // engine->Prealloc(max_length * 1.0, 5);
      return engine;
    };

    auto get_local_alignment_engine = [&](const std::size_t max_length) {
      auto engine = spoa::AlignmentEngine::Create(
          spoa::AlignmentType::kSW, // Smith-Waterman(local alignment)
          5,                        // match (default parameter form SPOA)
          -4,                       // mismatch
          -12,                      // gap
          -8                        // gap extension
      );
      // engine->Prealloc(max_length * 1.0, 5);
      return engine;
    };

    auto init_msa =
        [&](Window& w,
            std::shared_ptr<spoa::AlignmentEngine> alignment_engine) {
          assert(alignment_engine != nullptr);
          auto order = std::vector<std::size_t>(w.overlap_seqs.size());
          std::iota(order.begin(), order.end(), 0);
          std::random_shuffle(order.begin(), order.end());
          // alignment -> [(node_id_in_graph, id_in_seq), (), ...]
          // (id == -1) -> insertion or deletion
          auto alignment = alignment_engine->Align(w.backbone_seq, w.graph);
          w.graph.AddAlignment(alignment, w.backbone_seq, 1);
          for (const auto &idx : order) {
            alignment = alignment_engine->Align(w.overlap_seqs[idx], w.graph);
            /* default weight = 1 */
            w.graph.AddAlignment(alignment, w.overlap_seqs[idx], 1);
          }
          // ? The begin and end may need adjust due to sequencing error
          // ? happen at the begin and end of the read.
          // ? possible solution
          // ? 1. calcuate the coverage for first 10 and last 10 base, select
        };

    auto align_on_overlap_seqs =
        [&](Window& w, std::size_t read_idx,
            std::shared_ptr<spoa::AlignmentEngine>& alignment_engine) {
          assert(alignment_engine != nullptr);

          auto align_and_push_seq = [&](std::size_t idx, const std::string& seq,
                                        bool strain) {
            auto alignment = alignment_engine->Align(seq, w.graph);
            if (alignment.empty()) {
              return;
            }
            auto corrected = w.graph.DecodeAlignment(alignment);
            if (!strain) {
              corrected = bio::Codec::rev_comp(corrected);
            }
            std::scoped_lock lock(mutexes[idx]);
            reads[idx].corrected_fragments.emplace_back(std::move(corrected));
          };
          for (auto i = 0u; i < w.overlap_seqs.size(); i++) {
            align_and_push_seq(w.overlap_read_idx[i], w.overlap_seqs[i],
                               w.overlap_strains[i]);
          }
          align_and_push_seq(read_idx, w.backbone_seq, true);
        };

    auto batch_job = [&](Read& read) {
      thread_local std::shared_ptr<spoa::AlignmentEngine> local_alignment_engine =
          get_local_alignment_engine(Window::window_len);
      thread_local std::shared_ptr<spoa::AlignmentEngine> global_alignment_engine =
          get_global_alignment_engine(Window::window_len);
      
      auto build_window = [&](Window& w) {
        init_msa(w, global_alignment_engine);
        w.graph = w.graph.PruneGraph(w.get_prune_len());
        align_on_overlap_seqs(w, read.idx, local_alignment_engine);
      };
      for (auto& w : read.windows) {
        build_window(w);
      }
      {
        static std::atomic_int cnt = 0;
        spdlog::info("cnt = {}, read {} done, size = {}", ++cnt, read.idx,
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
    //       //                 w.idxL, std::chrono::duration_cast<std::chrono::milliseconds>(ed - st).count());
    //       //                 fs::path path = fmt::format(
    //       //       "/mnt/ec/ness/yolkee/thesis/tests/tmp/prune/{}", read.name);
    //       //   w.graph.PrintDot(path / fmt::format("{}-{}.dot", w.idxL, w.idxR));
    //       // }
    //     }
    //     static std::atomic_int cnt = 0;
    //     spdlog::info("cnt = {}, read {} done, size = {}", ++cnt, read.idx,
    //                  read.windows.size());
    //   }
    // };

#pragma omp parallel for num_threads(threads)
    for (auto& read : reads) {
      read.windows = make_windows(read);
    }
    std::vector<std::future<void>> futures;
    for (auto& read : reads) {
      auto [_, res] = threadpool.submit(batch_job, std::ref(read));
      futures.emplace_back(std::move(res));
    }

    // for (auto prev_idx = 0, idx = 0, sum = 0; idx < (int) reads.size(); idx++) {
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

    for (auto& f : futures) {
      f.get();
    }


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

    std::size_t fragment_size_sum = 0;
    #pragma omp parallel for reduction(+:fragment_size_sum)
    for (auto& read : reads) {
      std::size_t len_sum = 0;
      for (auto& fragment : read.corrected_fragments) {
        len_sum += fragment.size();
      }
      fragment_size_sum += read.corrected_fragments.size();
      spdlog::info("read {} corrected_fragments size = {}, average len = {}, "
                   "total len = {}, origin read len = {}",
                   read.idx, read.corrected_fragments.size(),
                   len_sum / read.corrected_fragments.size(), len_sum,
                   read.seq.size());
    }
    spdlog::info("average corrected_fragments size = {}",
                 fragment_size_sum / reads.size());

    auto path = fs::path("/mnt/ec/ness/yolkee/thesis/data/Ecoli/K12/seqs");
    fs::create_directories(path);
    for (auto& read : reads | std::views::take(10)) {
      auto p = path / (read.name + ".txt");
      spdlog::info("write to {}", p.string());
      std::ofstream fout(p);
      assert(fout.is_open() && "cannot open file");
      for (auto& fragment : read.corrected_fragments) {
        fout << fragment << '\n';
      }
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

    auto corrected_read = std::vector<bio::FastaRecord<false>>{};
    return corrected_read;
  }

  FragmentedReadCorrector(
      std::vector<bio::FastaRecord<false>>&& raw_reads,
      std::vector<bio::PafRecord>&& overlaps, const std::string& platform,
      const int thread_num = std::thread::hardware_concurrency(),
      bool debug = false)
      : BaseReadCorrector(std::move(raw_reads), std::move(overlaps), platform,
                          thread_num, debug),
        mutexes(raw_read_size) {
    make_reads();
    assemblers.reserve(raw_read_size);
    for (auto i : std::views::iota(0ul, raw_read_size)) {
      assemblers.emplace_back(ReadAssembler{ std::log2( reads[i].seq.size() )});
    }
  }

private:
  std::vector<ReadAssembler> assemblers;
  std::vector<Read> reads;
  std::vector<std::mutex> mutexes;
};