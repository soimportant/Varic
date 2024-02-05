#pragma once

#include <edlib.h>
#include <gperftools/heap-profiler.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <execution>
#include <future>
#include <mutex>
#include <ranges>
#include <spoa/spoa.hpp>
#include <string_view>
#include <thread>
#include <vector>

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

    int match = 5;
    int mismatch = -4;
    int gap = -8;
    int extend = -6;
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
  auto make_windows_for_one_read(Read& read) {
    read.init_windows(param.max_window_len, param.window_extend_len);
    for (const auto& overlap : read.overlap_range) {
      if (reads[overlap.t_id].need_corrected) {
        read.add_overlap_into_window(overlap, reads[overlap.t_id]);
      }
    }
    for (auto& t_id : read.overlap_reads_id) {
      auto fragments_cnt = read.get_fragment_cnt(t_id);
      auto index_range = reads[t_id].get_index_range(fragments_cnt);
      read.index_range[t_id] = index_range;
      read.indexes[t_id] = index_range.first;
    }
  }


  auto make_windows_for_all_read() {
    auto total_windows = 0ul;
    auto total_seqs = 0ul;
    auto total_seqs_lens = 0ul;
    const int take_reads = reads.size();
    // const int take_reads = 20;

    spdlog::debug("Taking {} reads now", take_reads);
    spdlog::info("Making windows for each reads...");

#pragma omp parallel for num_threads(threads) \
    reduction(+ : total_windows, total_seqs, total_seqs_lens)
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
    for (int i = 0; i < take_reads; i++) {
      reads[i].set_fragments_size();
    }

    spdlog::info("building window down, total = {}", total_windows);
    spdlog::debug("average windows in a read = {:.2f}",
                  static_cast<double>(total_windows) / take_reads);
    spdlog::debug("average sequence inside a window = {:.2f}",
                  static_cast<double>(total_seqs) / total_windows);
    spdlog::debug("average sequence length inside a window = {:.2f}",
                  static_cast<double>(total_seqs_lens) / total_seqs);
  }


 public:
  auto set_alignment_params(int match, int mismatch, int gap, int extend) {
    param.match = match;
    param.mismatch = mismatch;
    param.gap = gap;
    param.extend = extend;
  }

  /**
   * @brief correct the read
   *
   * @return auto
   */
  auto correct() {
    assert(std::ranges::is_sorted(overlaps, {}, &Overlap::q_id) &&
           "overlaps should be sorted by query id");

    auto cv = std::condition_variable{};

    // {
    //   auto file = fs::path("/mnt/ec/ness/yolkee/thesis/tmp/wrong.txt");
    //   auto fin = std::ifstream(file);
    //   auto read_names = std::vector<std::string>{};
    //   auto sampled_read_names = std::vector<std::string>{};
    //   auto name = std::string{};
    //   while (fin >> name) {
    //     read_names.emplace_back(std::move(name));
    //   }
    //   // random sample 150 reads
    //   std::ranges::sample(read_names, std::back_inserter(sampled_read_names),
    //                       150, std::mt19937{std::random_device{}()});
    //   for (auto& name : sampled_read_names) {
    //     auto read_id = name2id[name];
    //     wanted_read.emplace(read_id);
    //   }
    // }
    

    auto assemble_and_write = [&](Read& read) {
      // TODO: if coverage is enough, then we can assemble the read
      // remember that if we don't build the graph for this read, then there may
      // have other read that doesn't have enough data for building graph, so we
      // may need other mechanism to build the graph for this read.

      if (!read.ready_to_assemble()) {
        spdlog::debug("Read {} is not ready to assemble", read.name);
        return false;
      }
      // this check is not correct, the read.overlap_reads_id may not contain 
      // all source of fragments, due to preprocess step
      // for (auto t_id : read.overlap_reads_id) {
      //   if (!reads[t_id].ready_to_assemble() && !reads[t_id].is_assembled()) {
      //     spdlog::debug("Read {}'s overlap read {} is not ready to assemble",
      //                   read.name, reads[t_id].name);
      //     return false;
      //   }
      // }
      if (!read.acquire_assembling()) {
        spdlog::debug("Read {} acquire assembling failed", read.name);
        return false;
      }

      read.assemble_corrected_seq();
      static std::atomic_int cnt = 0;
      if (cnt % 100 == 0) {
        spdlog::debug("cnt = {}, all = {}", cnt, filtered_raw_reads_size);
      }
      cnt++;
      cv.notify_all();
      return true;
    };

    auto window_pipeline = [&](Window& window) {
      // TODO: check read is assembled or not, is yes, then skip this window

      auto read_id = window.read_id;
      // auto thread_id = threadpool.get_worker_id();

      thread_local auto global_aln_engine = get_global_alignment_engine(
        param.max_window_len,
        param.match,
        param.mismatch,
        param.gap,
        param.extend
      );
      thread_local auto semi_global_aln_engine = get_semi_global_alignment_engine(
        param.max_window_len,
        param.match,
        param.mismatch,
        param.gap,
        param.extend
      );
      auto corrected_sequences = window.get_corrected_sequences(global_aln_engine, semi_global_aln_engine);
      // window.build_variation_graph(engine);
      // {
      //   auto target_read = "01_10001";
      //   auto found = false;
      //   for (auto& seq : corrected_fragments) {
      //     if (seq.read_id == name2id[target_read]) {
      //       found = true;
      //     }
      //   }
      //   if (found) {
      //     auto path =
      //         fs::path(fmt::format("{}/corrected_fragments/{}_{}.txt",
      //         TMP_PATH,
      //                              id2name[window.read_id], window.idx));
      //     auto fout = std::ofstream(path);
      //     for (auto& seq : corrected_fragments) {
      //       auto fasta = bio::FastaRecord<false>{};
      //       fasta.name = id2name[seq.read_id];
      //       fasta.seq = seq.seq;
      //       fout << fasta << '\n';
      //     }
      //     window.print_graph(path.replace_extension(".dot"),
      //     name2id[target_read]);
      //   }
      // }
      for (auto& seq : corrected_sequences) {
        // auto lock = std::scoped_lock(mutexes[seq.read_id]);
        auto index = reads[read_id].get_index(seq.read_id);
        reads[seq.read_id].set_corrected_fragment(index, std::move(seq));
      }
      reads[read_id].add_finished_window();
      window.clear();
      // if (reads[window.read_id].ready_to_assemble()) {
      //   // all windows in this read has been processed, we can check there's
      //   // any read or itself need to be assembled
      //   threadpool.submit(assemble_and_write,
      //   std::ref(reads[window.read_id])); for (const auto& t_id :
      //   reads[window.read_id].overlap_reads_id) {
      //     threadpool.submit(assemble_and_write, std::ref(reads[t_id]));
      //   }
      // }
      {
        static int cnt = 0;
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
      // threadpool.submit(assemble_and_write, std::ref(read));

      if (read.id == name2id["01_10001"]) {
        spdlog::debug("read {} has {} windows", read.name, read.windows.size());
        for (auto i = 0u; i < read.windows.size(); i++) {
          spdlog::debug("window {}({}, {}) has {} overlap seqs", i,
                        read.windows[i].start, read.windows[i].end,
                        read.windows[i].overlap_seqs.size());
        }
      }
    };



    // spdlog::debug("Pushed all read into pipeline...");
    // // TODO fixed the order
    // auto order = std::vector<std::size_t>(reads.size());
    // std::iota(order.begin(), order.end(), 0);
    // std::sort(order.begin(), order.end(), [&](int a, int b) {
    //   return reads[a].overlap_range.size() < reads[b].overlap_range.size();
    // });

    // for (auto& x : order) {
    //   auto& read = reads[x];
    //   if (read.need_corrected) {
    //     auto [_, res] = threadpool.submit(read_pipeline, std::ref(read));
    //   }
    // }
    // spdlog::debug("Pushed read pipeline done");

    {
      make_windows_for_all_read();
      spdlog::info("Build windows for all reads");

      auto futures = std::vector<std::future<void>>{};
      for (auto i = 0u; i < reads.size(); i++) {
        auto& read = reads[i];
        // {
        //   bool found = false;
        //   if (wanted_read.contains(read.id)) {
        //     found = true;
        //   }
        //   for (auto& t_id : read.overlap_reads_id) {
        //     if (wanted_read.contains(t_id)) {
        //       found = true;
        //     }
        //   }
        //   if (!found) {
        //     continue;
        //   }
        // }
        for (auto& w : read.windows) {
          auto [_, res] = threadpool.submit(window_pipeline, std::ref(w));
          futures.emplace_back(std::move(res));
        }
      }

      for (auto& res : futures) {
        res.get();
      }
      auto total_windows = 0ul;
      for (auto& read : reads) {
        total_windows += read.windows.size();
      }
      spdlog::info("Build windows done: total = {}", total_windows);
    }

    {
      // spdlog::debug("check all window is ready for assembled");
      // for (auto& read : reads) {
      //   if (read.need_corrected && !read.ready_to_assemble()) {
      //     spdlog::debug("read {} is not ready to assemble", read.name);
      //   }
      //   for (auto& t_id : read.overlap_reads_id) {
      //     if (!reads[t_id].ready_to_assemble()) {
      //       spdlog::debug("Read {} is not ready to assemble: {}",
      //                     reads[t_id].name, reads[t_id].need_corrected);
      //     }
      //   }
      // }
    }

    {
      spdlog::info("Assembling all reads");
      auto futures = std::vector<std::future<bool>>{};
      for (auto i = 0u; i < reads.size(); i++) {
        auto& read = reads[i];
        // {
        //   bool found = false;
        //   if (wanted_read.contains(read.id)) {
        //     found = true;
        //   }
        //   if (!found) {
        //     continue;
        //   }
        // }
        if (read.need_corrected && read.ready_to_assemble()) {
          auto [_, res] = threadpool.submit(assemble_and_write, std::ref(read));
          futures.emplace_back(std::move(res));
        }
      }
      for (auto& res : futures) {
        assert(res.get());
      }
    }

    // threadpool.join();
    // spdlog::debug("After join, corrected_reads.size() = {}",
    // corrected_reads.size());

    // {
    //   spdlog::debug("Waiting all reads been assembled");
    //   auto mutex = std::mutex{};
    //   std::unique_lock lock(mutex);
    //   cv.wait(lock, [&]() {
    //     spdlog::debug("cv got notified, corrected_reads.size() = {}",
    //     corrected_reads.size()); return corrected_reads.size() ==
    //     filtered_raw_reads_size;
    //   });
    //   spdlog::debug("Done");
    // }

    // Do this check for there are some reads that are not assembled

    // std::exit(0);

    std::vector<bio::FastaRecord<false>> corrected_reads;

    // TODO: remove this, use below
    // for (auto id : wanted_read) {
    //   auto& read = reads[id];
    //   if (read.need_corrected) {
    //     auto corrected_read = bio::FastaRecord<false>{
    //         .name = read.name,
    //         .seq = std::move(read.corrected_seq),
    //     };
    //     if (corrected_read.seq.size() != 0) {
    //       corrected_reads.emplace_back(std::move(corrected_read));
    //     }
    //   }
    // }

#pragma omp parallel for num_threads(threads)
    for (auto& read : reads) {
      // spdlog::debug("Assemble read {}, len = {}", id2name[read.id],
      // read.len());
      auto corrected_read = bio::FastaRecord<false>{
          .name = read.name,
          .seq = std::move(read.corrected_seq),
      };
#pragma omp critical
      if (corrected_read.seq.size() != 0) {
        corrected_reads.emplace_back(std::move(corrected_read));
      }
    }
    spdlog::debug("Filtered reads = {}", reads.size());
    spdlog::debug("corrected_reads.size() = {}", corrected_reads.size());
    return corrected_reads;
  }

private:
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
      return true;
    };

    reads.resize(raw_reads.size());
#pragma omp parallel for
    for (auto& r : raw_reads) {
      auto read = Read(std::move(r));
      read.id = name2id[read.name];
      read.need_corrected = valid_read(read);
      reads[read.id] = std::move(read);
    }

    auto need_create_rc = std::vector<std::atomic_bool>(reads.size());
#pragma omp parallel for
    for (auto& read : reads) {
      read.overlap_range = get_overlap_range(read.name);
      for (const auto& overlap : read.overlap_range) {
        if (reads[overlap.t_id].need_corrected) {
          read.overlap_reads_id.emplace_back(overlap.t_id);
        }
        if (!overlap.forward_strain) {
          need_create_rc[read.id] = true;
          need_create_rc[overlap.t_id] = true;
        }
      }
    }

    for (auto i = 0u; i < reads.size(); i++) {
      if (need_create_rc[i]) {
        reads[i].create_rc();
      }
    }

    for (auto& read : reads) {
      if (!read.check()) {
        spdlog::error("read {} check failed", read.name);
        exit(EXIT_FAILURE);
      }
      for (auto& overlap : read.overlap_range) {
        if (!overlap.forward_strain && !reads[overlap.t_id].check()) {
          spdlog::error("overlap {} <-> {} check failed", read.name,
                        reads[overlap.t_id].name);
        }
      }
    }
    spdlog::debug("Pass precheck");

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
      auto valid_read = [&](const std::string& name) {
        return name2id.contains(name);
      };
      if (!valid_read(paf.q_name) || !valid_read(paf.t_name) ||
          paf.q_name == paf.t_name) {
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

  auto check() {
    // auto found = [&](const int& q_id, const int& t_id) {
    //   for (auto& overlap : reads[q_id].overlap_range) {
    //     if (overlap.t_id == t_id) {
    //       return true;
    //     }
    //   }
    //   return false;
    // };
    // for (auto& read : reads) {
    //   if (!wanted_read.contains(read.id)) {
    //     continue;
    //   }
    //   for (auto t_id : read.overlap_reads_id) {
    //     if (!found(t_id, read.id)) {
    //       spdlog::error("read {} <-> {} not found", read.name,
    //                     reads[t_id].name);
    //     }
    //   }
    // }
  }

  auto preprocess(std::vector<R>& raw_reads,
                  std::vector<bio::PafRecord>& overlaps) {
    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();
    mutexes = std::vector<std::mutex>(raw_reads.size());
    /* transform read name to id */
    for (auto& read : raw_reads) {
      auto id = name2id.size();
      name2id[read.name] = id;
      id2name[id] = read.name;
    }
    overlap_preprocess(overlaps);
    read_preprocess(raw_reads);
    check();
  }

public:
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
      : threads(thread_num),
        threadpool(bio::make_threadpool(this->threads)),
        platform(platform),
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
  std::mutex corrected_reads_mutex;

  /* debug flag */
  bool debug = false;

  std::set<int> wanted_read;
};