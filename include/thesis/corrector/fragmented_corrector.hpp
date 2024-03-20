#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdlib>
#include <execution>
#include <future>
#include <mutex>
#include <random>
#include <ranges>
#include <string_view>
#include <thread>
#include <vector>

#include <edlib.h>
#include <gperftools/heap-profiler.h>
#include <omp.h>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/dynamic_progress.hpp>
#include <indicators/multi_progress.hpp>
#include <indicators/progress_bar.hpp>
#include <spoa/spoa.hpp>

#include "indicators/indeterminate_progress_bar.hpp"
#include "thesis/corrector/detail/overlap.hpp"
#include "thesis/corrector/detail/read.hpp"
#include "thesis/corrector/detail/window.hpp"
#include "thesis/format/paf.hpp"

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
  using id_t = typename Read::id_t;

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

    /* maximum used sequence when building variation graph */
    int depth = -1;

    double pruned_ratio = 0.95;

    std::int64_t seed = 0;
  } param;

  /**
   * @brief Initializes the progress bars for different stages of the correction
   * process.
   *
   * This function sets up progress bars using the indicators library for
   * tracking the progress of various stages in the correction process. It
   * configures the appearance and options of each progress bar, such as bar
   * width, fill character, prefix text, color, and font style. The progress
   * bars are used to display the progress of preprocessing, window
   * initialization, and collecting corrected sequences.
   */
  auto init_progress_bar() {
    using namespace indicators;
    auto setup = [&](auto& bar, std::string prefix, auto color = Color::white) {
      prefix = fmt::format("{:<30}", prefix);
      bar.set_option(option::BarWidth{30});
      bar.set_option(option::Start{" ["});
      bar.set_option(option::Fill{"="});
      bar.set_option(option::Lead{">"});
      bar.set_option(option::Remainder{" "});
      bar.set_option(option::End{" ] "});
      bar.set_option(option::PrefixText{prefix});
      bar.set_option(option::ForegroundColor{color});
      bar.set_option(option::ShowElapsedTime{true});
      bar.set_option(option::ShowRemainingTime{true});
      bar.set_option(option::ShowPercentage(true));
      bar.set_option(
          option::FontStyles{std::vector<FontStyle>{FontStyle::bold}});
    };

    setup(preprocess_bar, "Preprocess", Color::yellow);
    setup(make_window_bar, "Init windows", Color::yellow);
    setup(collect_corrected_seq_bar, "Collect corrected sequences",
          Color::yellow);
    setup(assemble_corrected_seq_bar, "Assembled corrected read",
          Color::yellow);
    bars.push_back(collect_corrected_seq_bar);
    bars.push_back(assemble_corrected_seq_bar);
  }

  /**
   * @brief Creates windows for a single read.
   *
   * This function initializes windows for a given read based on the maximum
   * window length and window extend length parameters. It then adds overlaps
   * from other reads into the windows if the corresponding reads need to be
   * corrected. Finally, it calculates the index range and indexes for each
   * overlap read.
   *
   * @param read The read for which windows are created.
   */
  auto make_windows_for_one_read(Read& read) {
    read.init_windows(param.max_window_len, param.window_extend_len,
                      param.seed);
    for (const auto& overlap : read.overlap_range) {
      if (reads[overlap.t_id].need_corrected) {
        read.add_overlap_into_window(overlap, reads[overlap.t_id]);
      }
    }
    for (auto& w : read.windows) {
      w.shrink_to_fit();
    }
    for (auto& t_id : read.overlap_reads_id) {
      auto fragments_cnt = read.get_fragment_cnt(t_id);
      auto index_range = reads[t_id].get_index_range(fragments_cnt);
      read.index_range[t_id] = index_range;
      read.indexes[t_id] = index_range.first;
    }
  }

  auto make_windows_for_all_read() {
    std::size_t total_windows = 0ul;
    std::size_t total_seqs = 0ul;
    std::size_t total_seqs_lens = 0ul;
    const int take_reads = reads.size();
    // const int take_reads = 20;

    spdlog::debug("Taking {} reads now", take_reads);
    indicators::show_console_cursor(false);
    make_window_bar.set_option(indicators::option::MaxProgress(take_reads));

#pragma omp parallel num_threads(threads)
    {
#pragma omp for reduction(+ : total_windows, total_seqs, total_seqs_lens)
      for (int i = 0; i < take_reads; i++) {
        auto& read = reads[i];
        make_windows_for_one_read(read);
        total_windows += read.windows.size();
        for (auto& w : read.windows) {
          total_seqs += w.overlap_seqs.size();
          for (auto& seq : w.overlap_seqs) {
            total_seqs_lens += seq.seq.size();
          }
        }
        make_window_bar.tick();
      }
      make_window_bar.set_option(indicators::option::PostfixText{
          "allocating space for corrected sequences"});
#pragma omp for
      for (int i = 0; i < take_reads; i++) {
        reads[i].set_fragments_size();
      }
    }
    collect_corrected_seq_bar.set_option(
        indicators::option::MaxProgress(total_windows));

    spdlog::info("building window down, total = {}", total_windows);
    spdlog::debug("average windows in a read = {:.2f}",
                  static_cast<double>(total_windows) / take_reads);
    spdlog::debug("average sequence inside a window = {:.2f}",
                  static_cast<double>(total_seqs) / total_windows);
    spdlog::debug("average sequence length inside a window = {:.2f}",
                  static_cast<double>(total_seqs_lens) / total_seqs);
    auto dummy = std::vector<Overlap>{};
    std::swap(overlaps, dummy);
    return total_windows;
  }

  /**
   * @brief Performs read preprocessing.
   *
   * This function processes a vector of raw reads by performing various checks
   * and modifications on each read.
   *
   * @param raw_reads The vector of raw reads to be preprocessed.
   */
  auto read_preprocess(std::vector<R> raw_reads) {
    preprocess_bar.set_option(
        indicators::option::PostfixText{"Read preprocessing"});

    reads.resize(raw_reads.size());
    auto need_create_rc = std::vector<std::atomic_bool>(reads.size());

    auto identify_valid_read = [&]() {
      auto valid_read = [&](Read& r) {
        if (r.seq.size() < param.min_read_length) {
          return false;
        }
        return true;
      };
#pragma omp parallel for num_threads(threads)
      for (auto& r : raw_reads) {
        auto read = Read(std::move(r));
        read.id = name2id[read.name];
        read.need_corrected = valid_read(read);
        reads[read.id] = std::move(read);
      }
      filtered_raw_reads_size =
          std::ranges::count_if(reads, &Read::need_corrected);
      assemble_corrected_seq_bar.set_option(
          indicators::option::MaxProgress(filtered_raw_reads_size));
    };

    auto identify_overlap_range = [&]() {
      auto get_overlap_range = [&](const std::string& name) {
        auto id = name2id[name];
        auto st = std::ranges::lower_bound(overlaps, id, {}, &Overlap::q_id);
        auto ed = std::ranges::upper_bound(overlaps, id, {}, &Overlap::q_id);
        return std::ranges::subrange(st, ed);
      };
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
    };

    auto create_rc = [&] {
      for (auto i = 0u; i < reads.size(); i++) {
        if (need_create_rc[i]) {
          reads[i].create_rc();
        }
      }
    };

    identify_valid_read();
    preprocess_bar.set_progress(80);
    identify_overlap_range();
    preprocess_bar.set_progress(95);
    create_rc();
    preprocess_bar.set_progress(100);

    // for (auto& read : reads) {
    //   if (!read->check()) {
    //     spdlog::error("read {} check failed", read->name);
    //     exit(EXIT_FAILURE);
    //   }
    //   for (auto& overlap : read->overlap_range) {
    //     if (!overlap.forward_strain && !reads[overlap.t_id].check()) {
    //       spdlog::error("overlap {} <-> {} check failed", read->name,
    //                     reads[overlap.t_id]->name);
    //     }
    //   }
    // }
  }

  /**
   * @brief
   * @param raw_reads
   * @param overlaps
   * @return
   */
  auto overlap_preprocess(std::vector<bio::PafRecord>&& overlaps) {
    preprocess_bar.set_option(
        indicators::option::PostfixText{"Overlap preprocessing"});

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
    filtered_overlaps.reserve(overlaps.size());
#pragma omp parallel for num_threads(threads)
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
    preprocess_bar.set_progress(70);
  }

  /**
   * @brief debug function
   *
   * @return auto
   */
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

  /**
   * Preprocesses the raw reads and overlaps.
   *
   * @param raw_reads The vector of raw reads.
   * @param overlaps The vector of overlaps.
   */
  auto preprocess(std::vector<R>&& raw_reads,
                  std::vector<bio::PafRecord>&& overlaps) {
    init_progress_bar();
    unfiltered_raw_read_size = raw_reads.size();
    unfiltered_overlap_size = overlaps.size();

    indicators::show_console_cursor(false);
    preprocess_bar.set_progress(0);

    /* transform read name to id */
    for (auto& read : raw_reads) {
      auto id = name2id.size();
      name2id[read.name] = id;
      id2name[id] = read.name;
    }

    preprocess_bar.set_progress(5);
    preprocess_bar.set_option(
        indicators::option::PostfixText{"Transforming read name to id"});
    overlap_preprocess(std::move(overlaps));
    read_preprocess(std::move(raw_reads));
    check();
  }

  /**
   * @brief Collects corrected fragments for each read.
   *
   * This function collects corrected fragments for each read by processing
   * windows in parallel. It uses global and local alignment engines to correct
   * the fragments and updates the reads with the corrected fragments. Once all
   * windows have been processed for a read, the function checks if the read is
   * ready to be assembled. If so, it submits the read for assembly and writing.
   * The progress of the collection process is tracked using a progress bar.
   *
   * @note This function assumes that the necessary data structures and
   * parameters have been initialized before calling it.
   */
  auto collect_corrected_fragments() {
    auto cv = std::condition_variable{};
    static std::atomic_int finished_windows_cnt = 0;
    static std::atomic_int finished_reads_cnt = 0;

    auto window_pipeline = [&](Window& window) {
      auto& read = reads[window.read_id];

      thread_local auto global_aln_engine =
          get_global_alignment_engine(param.max_window_len, param.match,
                                      param.mismatch, param.gap, param.extend);
      thread_local auto local_aln_engine =
          get_local_alignment_engine(param.max_window_len, param.match,
                                     param.mismatch, param.gap, param.extend);

      // mask indicate which read should corrected sequence
      auto read_mask = std::vector<bool>(window.depth(), true);
      for (auto i = 0; i < window.depth(); i++) {
        auto t_id = window.overlap_seqs[i].read_id;
        if (!reads[t_id].need_corrected || reads[t_id].is_assembled()) {
          read_mask[i] = false;
        }
      }
      auto corrected_fragments = window.get_corrected_fragments(
          global_aln_engine, local_aln_engine, read_mask, param.pruned_ratio,
          param.depth);
      for (auto& seq : corrected_fragments) {
        auto& overlap_read = reads[seq.read_id];
        if (overlap_read.is_assembled()) {
          continue;
        }
        auto index = read.get_index(seq.read_id);
        overlap_read.set_corrected_fragment(index, std::move(seq));
        if (overlap_read.is_assembled()) {
          bars[1].tick();
          bars[1].set_option(indicators::option::PostfixText{
              fmt::format("{:2d}/{:6d}", bars[1].current(),
              filtered_raw_reads_size)});
        }
      }
      read.add_finished_window();
      finished_windows_cnt++;
      if (finished_windows_cnt % 50 == 0) {
        bars[0].set_progress(finished_windows_cnt);
        // collect_corrected_seq_bar.set_progress(finished_windows_cnt);
      }
      cv.notify_one();
      // if (reads[window.read_id].ready_to_assemble()) {
      //   // all windows in this read has been processed, jwe can check there's
      //   // any read or itself need to be assembled
      //   threadpool.submit(assemble_and_write,
      //   std::ref(reads[window.read_id]));
      // for (const auto& t_id :
      //   reads[window.read_id].overlap_reads_id) {
      //     threadpool.submit(assemble_and_write, std::ref(reads[t_id]));
      //   }
      // }
    };
    auto total_windows = make_windows_for_all_read();
    // std::exit(0);

    // spdlog::info("Build windows for all reads");
    for (auto i = 0u; i < reads.size(); i++) {
      auto& read = reads[i];
      for (auto& w : read.windows) {
        boost::asio::post(
            threadpool, [&, w = std::ref(w)]() mutable { window_pipeline(w); });
      }
    }
    std::mutex mutex;
    std::unique_lock lock(mutex);
    cv.wait(lock, [&]() {
      return finished_windows_cnt == total_windows;
      // return collect_corrected_seq_bar.is_completed();
    });
    bars[0].mark_as_completed();
  }

  auto assemble_corrected_read() {
    auto assemble_and_write = [&](Read& read) {
      // this check is not correct, the read.overlap_reads_id may not contain
      // all source of fragments, due to preprocess step
      // for (auto t_id : read.overlap_reads_id) {
      //   if (!reads[t_id].ready_to_assemble() && !reads[t_id].is_assembled())
      //   {
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
      // add_corrected_read_and_reset(read.id);
      bars[1].tick();
      bars[1].set_option(indicators::option::PostfixText{
          fmt::format("{:2d}/{:6d}", bars[1].current(),
                      filtered_raw_reads_size)});
      return true;
    };

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
      if (read.need_corrected && read.is_assembled()) {
        continue;
      }
      boost::asio::post(threadpool, [&, read = std::ref(read)]() mutable {
        assemble_and_write(read);
      });
    }
    threadpool.wait();
  }

 public:
  /**
   * @brief Sets the alignment parameters for the corrector.
   *
   * @param match The score for a match between two characters.
   * @param mismatch The score for a mismatch between two characters.
   * @param gap The score for introducing a gap in the alignment.
   * @param extend The score for extending an existing gap in the alignment.
   */
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
    collect_corrected_fragments();
    // {
    //   auto dir = fs::path("/mnt/ec/ness/yolkee/thesis/tmp/fragments");
    //   if (!fs::exists(dir)) {
    //     fs::create_directories(dir);
    //   } else {
    //     fs::remove_all(dir);
    //     fs::create_directories(dir);
    //   }
    //   for (auto i = 0u; i < reads.size(); i++) {
    //     auto& read = reads[i];
    //     if (read.need_corrected) {
    //       auto path = dir / fmt::format("{}.txt", read.name);
    //       read.save_corrected_fragments(path);
    //     }
    //   }
    // }
    assemble_corrected_read();

    std::vector<bio::FastaRecord<false>> corrected_reads;

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
    spdlog::info("Thread number: {}", threads);
    spdlog::info("(Match, Mismatch, Gap, Gap Extend) = ({}, {}, {}, {})",
                 param.match, param.mismatch, param.gap, param.extend);
    spdlog::info("Used sequence when build variation graph: {}", param.depth);
    spdlog::info("Variation graph pruning ratio: {:.3f}", param.pruned_ratio);
    spdlog::info("Random seed: {}", param.seed);
  }

  FragmentedReadCorrector(
      std::vector<R>&& raw_reads, std::vector<bio::PafRecord>&& overlaps,
      const std::string& platform,
      const int thread_num = std::thread::hardware_concurrency(),
      const int depth = -1, const double pruned_ratio = 0.95,
      const std::int64_t seed = 0)
      : threads(thread_num), threadpool(this->threads), platform(platform) {
    preprocess(std::move(raw_reads), std::move(overlaps));
    param.depth = depth;
    param.pruned_ratio = pruned_ratio;
    param.seed = seed;
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
  // std::vector<std::unique_ptr<Read>> reads;
  std::vector<Read> reads;
  std::vector<std::mutex> mutexes;

  /* overlap information between `raw_reads` */
  std::vector<Overlap> overlaps;

  /* mapping from read name to read id and vice versa */
  std::map<std::string, id_t, std::less<>> name2id;
  std::map<id_t, std::string> id2name;

  /* max threads and threadpool */
  std::size_t threads;
  boost::asio::thread_pool threadpool;

  /* sequencing platform of raw_reads */
  // TODO: should be enum class
  std::string platform;
  std::mutex corrected_reads_mutex;
  std::vector<bio::FastaRecord<false>> corrected_reads;

  /* random engine */
  std::mt19937_64 random_engine;

  /* progress bars */
  indicators::ProgressBar preprocess_bar;
  indicators::ProgressBar make_window_bar;
  indicators::ProgressBar collect_corrected_seq_bar;
  indicators::ProgressBar assemble_corrected_seq_bar;
  indicators::DynamicProgress<indicators::ProgressBar> bars;

  std::set<int> wanted_read;
};