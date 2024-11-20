#pragma once

#include <atomic>
#include <concepts>
#include <optional>
#include <queue>
#include <ranges>

#include <edlib.h>
#include <spdlog/spdlog.h>
#include <biovoltron/file_io/fasta.hpp>

#include "thesis/algo/assemble/read_assembler.hpp"
#include "thesis/corrector/detail/overlap.hpp"
#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/corrector/detail/window.hpp"
#include "thesis/ds/segtree.hpp"

namespace bio = biovoltron;

/**
 * @brief align overlap part between query read and target read, and return
 * the cigar string
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

template <class R>
  requires std::derived_from<R, bio::FastaRecord<R::encoded>>
class TestReadWrapper : public R {
 public:
  /* disable copy constructor and copy assignment */
  TestReadWrapper(const TestReadWrapper&) = delete;
  TestReadWrapper& operator=(const TestReadWrapper&) = delete;

  /* move constructor */
  TestReadWrapper(TestReadWrapper&& rhs)
      : R(std::move(rhs)), coverage(std::move(rhs.coverage)) {
    this->id = rhs.id;
    this->overlap_reads_id = std::move(rhs.overlap_reads_id);
    this->overlap_range = std::move(rhs.overlap_range);
    this->windows = std::move(rhs.windows);
    this->finished_windows_cnt = rhs.finished_windows_cnt.load();
    this->fragment_cnt = std::move(rhs.fragment_cnt);
    this->index_range = std::move(rhs.index_range);
    this->indexes = std::move(rhs.indexes);
    this->fragments = std::move(rhs.fragments);
    this->total_fragments = rhs.total_fragments.load();
    this->corrected_seq = std::move(rhs.corrected_seq);
    this->rc_seq = std::move(rhs.rc_seq);
    this->rev_qual = std::move(rhs.rev_qual);
    this->state = rhs.state.load();
  }

  /* move assignment */
  TestReadWrapper& operator=(TestReadWrapper&& rhs) {
    this->id = rhs.id;
    this->name = std::move(rhs.name);
    this->seq = std::move(rhs.seq);
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      this->qual = std::move(rhs.qual);
    }
    this->overlap_reads_id = std::move(rhs.overlap_reads_id);
    this->overlap_range = std::move(rhs.overlap_range);
    this->windows = std::move(rhs.windows);
    this->finished_windows_cnt = rhs.finished_windows_cnt.load();
    this->fragment_cnt = std::move(rhs.fragment_cnt);
    this->index_range = std::move(rhs.index_range);
    this->indexes = std::move(rhs.indexes);
    this->fragments = std::move(rhs.fragments);
    this->total_fragments = rhs.total_fragments.load();
    this->coverage = std::move(rhs.coverage);
    this->corrected_seq = std::move(rhs.corrected_seq);
    this->rc_seq = std::move(rhs.rc_seq);
    this->rev_qual = std::move(rhs.rev_qual);
    this->state = rhs.state.load();
    return *this;
  }

  TestReadWrapper(R&& r) : R(std::move(r)), coverage(this->seq.size()) {}

  TestReadWrapper() {}

  enum struct State {
    UNCORRECTED,
    BUILDING_WINDOW,
    WAIT_FOR_ASSEMBLE,
    ASSEMBLING,
    ASSEMBLED_SUCCESS,
    ASSEMBLED_FAILED,
  };

  /**
   * @brief Initialize boundaries and backbone sequence for each window.
   * @details The window length will be equally divided by the length of read,
   * and the max length of window is `max_window_len`. Then extend the window by
   * `window_extend_len` on both side.
   *
   * @param max_window_len max length of window
   * @param window_extend_len extend length of window
   * @return void
   */
  auto init_windows(const std::size_t max_window_len,
                    const std::size_t window_extend_len) {
    const auto window_sz =
        len() / max_window_len + (len() % max_window_len != 0);
    const auto window_len = len() / window_sz + (len() % window_sz != 0);

    for (auto pos = 0ul, w_idx = 0ul; pos < len(); pos += window_len, w_idx++) {
      auto window =
          Window{this->id, w_idx, pos, std::min(len(), pos + window_len)};
      window.extend(window_extend_len, len());

      window.backbone.read_id = this->id;
      window.backbone.left_bound = window.start;
      window.backbone.right_bound = window.end;
      std::tie(window.backbone.seq, window.backbone.qual) =
          subview(true, window.start, window.end);
      window.backbone.forward_strain = true;

      windows.emplace_back(std::move(window));
    }
    fragment_cnt[this->id] = windows.size();
    index_range[this->id] = get_index_range(windows.size());
    indexes[this->id] = index_range[this->id].first;
    total_fragments += windows.size();
    state = State::BUILDING_WINDOW;
  }

  /**
   * Retrieves a subrange of windows within a specified range.
   *
   * This function searches for the window range that contains the specified
   * start and end positions. It returns a subrange of windows that fall within
   * the specified range.
   *
   * @param start The start position of the range.
   * @param end The end position of the range.
   * @return A subrange of windows within the specified range.
   */
  auto get_window_range(std::size_t start, std::size_t end) {
    auto st = std::ranges::upper_bound(windows, start, {}, &Window::start);
    if (st != windows.end()) {
      st = std::prev(st);
    }
    auto ed = std::ranges::lower_bound(windows, end, {}, &Window::start);
    return std::ranges::subrange(st, ed);
  }

  /**
   * @brief Find breakpoints from a cigar string within a window range.
   *
   * This function takes a cigar string, a window range, and an overlap as
   * input. It iterates through the cigar string and determines the breakpoints
   * in the target read and query read for each window in the window range.
   *
   * @link https://github.com/isovic/racon/blob/master/src/overlap.cpp
   *
   * @param cigar The cigar string.
   * @param window_range The range of windows.
   * @param overlap The overlap between the target read and query read.
   * @return A pair of vectors containing the breakpoints in the query read and
   *         target read for each window.
   */
  template <std::ranges::range T>
  auto find_breakpoints_from_cigar(const bio::Cigar& cigar,
                                   const T window_range,
                                   const Overlap& overlap) {
    auto window_cnt = std::ranges::size(window_range);
    auto [q_start, q_end] = std::make_pair(overlap.q_idx_L, overlap.q_idx_R);
    auto [t_start, t_end] = std::make_pair(overlap.t_idx_L, overlap.t_idx_R);
    auto forward_strain = overlap.forward_strain;

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
     * @brief set the start position of target read when meet a match or
     * mismatch in cigar string for each window.
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

    for (auto w_idx = 0u; const auto& [len, op] : cigar) {
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
    if (!forward_strain) {
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

  /**
   * @brief Adds an overlap into the window.
   *
   * This function takes an overlap and a read wrapper as input and adds the
   * corresponding sequence into the window. It takes the `std::string_view` on
   * both query and target read and aligns them using `align` function. Then it
   * finds the breakpoints by the cigar string using
   * `find_breakpoints_from_cigar` function. The function ensures that the size
   * of q_breakpoints and t_breakpoints is equal, then iterates over the windows
   * and adds the sequence into each window according to the breakpoints.
   *
   * @param overlap The overlap to be added into the window.
   * @param t_read The read wrapper containing the sequence info to be added.
   */
  auto add_overlap_into_window(const Overlap& overlap,
                               const TestReadWrapper<R>& t_read) {
    const auto forward_strain = overlap.forward_strain;

    const auto [q_start, q_end] =
        std::make_pair(overlap.q_idx_L, overlap.q_idx_R);
    const auto [t_start, t_end] =
        std::make_pair(overlap.t_idx_L, overlap.t_idx_R);
    const auto [q_seq_view, q_qual_view] = subview(true, q_start, q_end);
    const auto [t_seq_view, t_qual_view] =
        t_read.subview(forward_strain, t_start, t_end);

    const auto window_range = get_window_range(q_start, q_end);
    auto window_idx = window_range.begin() - windows.begin();

    auto cigar = align(q_seq_view, t_seq_view);
    auto [q_breakpoints, t_breakpoints] =
        find_breakpoints_from_cigar(cigar, window_range, overlap);
    assert(q_breakpoints.size() == t_breakpoints.size() &&
           "size not equal between q_breakpoints and t_breakpoints");

    const auto window_cnt = std::ranges::size(window_range);
    for (auto i = 0u; i < window_cnt; i++, window_idx++) {
      auto& window = windows[window_idx];
      const auto [t_window_st, t_window_ed] = t_breakpoints[i];
      if (t_window_st >= t_end) {
        continue;
      }
      assert(t_window_ed >= t_window_st &&
             "t_window_ed < t_window_st, wrong boundary");
      const auto [t_window_seq, t_window_qual] =
          t_read.subview(forward_strain, t_window_st, t_window_ed);

      // when racon add sequence into window, it will store the pos of first
      // match, then sort the overlap_seqs according the pos of first match
      window.add_sequence(
          Sequence<>{
              overlap.t_id,   // read_id
              t_window_st,    // left_bound
              t_window_ed,    // right_bound
              t_window_seq,   // seq
              t_window_qual,  // qual
              forward_strain  // forward_strain
          },
          q_breakpoints[i]);
    }
    fragment_cnt[t_read.id] += window_cnt;
  }

  /**
   * @brief Increases the count of finished windows.
   *
   * @details Each time this function is called, the count of finished windows
   * will be increased by one. If the count of finished windows equals the total
   * number of windows, the status will be set to WAIT_TO_ASSEMBLE.
   */
  auto add_finished_window() noexcept {
    finished_windows_cnt++;
    if (finished_windows_cnt == windows.size()) {
      auto expected = State::BUILDING_WINDOW;
      state.compare_exchange_strong(expected, State::WAIT_FOR_ASSEMBLE);
      return true;
    }
    return false;
  }

  /**
   * @brief Returns whether the read is ready to assemble.
   *
   * @return true if the read is ready to assemble, false otherwise.
   */
  auto ready_to_assemble() noexcept {
    return state == State::WAIT_FOR_ASSEMBLE;
  }

  /**
   * Acquires the assembling state.
   * This function atomically compares the current state with
   * State::WAIT_FOR_ASSEMBLE, and if they are equal, it updates the state to
   * State::ASSEMBLING and returns true. If the current state is not equal to
   * State::WAIT_FOR_ASSEMBLE, it returns false.
   *
   * @return true if the state was successfully updated to State::ASSEMBLING,
   * false otherwise.
   */
  auto acquire_assembling() {
    auto expected = State::WAIT_FOR_ASSEMBLE;
    return state.compare_exchange_strong(expected, State::ASSEMBLING);
  }

  auto is_assembled() const noexcept { return state == State::ASSEMBLED_SUCCESS; }

  /**
   * @brief Creates the reverse complement sequence and quality scores (if
   * applicable).
   *
   * @return true if the reverse complement sequence is successfully created,
   * false otherwise.
   */
  auto create_rc() noexcept {
    if (rc_seq.size() != 0) {
      return;
    }
    rc_seq = bio::Codec::rev_comp(this->seq);
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      rev_qual = std::string(this->qual.rbegin(), this->qual.rend());
    }
  }

  /**
   * @brief Returns the length of the sequence.
   *
   * @return The length of the sequence.
   */
  auto len() const noexcept { return this->seq.size(); }

  /**
   * @brief Returns a subview of the sequence and quality data.
   *
   * This function returns a subview of the sequence and quality data based on
   * the specified start and end positions. The subview is determined by the
   * `forward_strain` parameter, which indicates whether to use the forward or
   * reverse complement sequence. The start and end positions define the range
   * of the subview.
   *
   * @param forward_strain A boolean value indicating whether to use the forward
   * strain or reverse complement sequence.
   * @param start The start position of the subview.
   * @param end The end position of the subview.
   * @return A pair containing the subview of the sequence data and an optional
   * subview of the quality data.
   * @throws std::invalid_argument If the start position is larger than the end
   * position.
   * @throws std::out_of_range If the start position is out of range of the
   * sequence.
   */
  auto subview(bool forward_strain, const std::size_t start,
               const std::size_t end) const {
    if (start > end) {
      spdlog::error("start position {} is behind end position", start, end);
      throw std::invalid_argument("start position is behind end position");
    }

    if (start > this->seq.size()) {
      spdlog::error(
          "start position {} is out of range, length = {}, strand = {}", start,
          this->seq.size(), forward_strain);
      throw std::out_of_range("start position is out of range of sequence");
    }

    auto real_start = forward_strain ? start : this->seq.size() - end;
    auto len = end - start;

    auto seq_view = std::string_view{};
    seq_view = forward_strain ? this->seq : rc_seq;
    if (seq_view.size() < real_start + len) {
      spdlog::error("forward strain = {}, real_start = {}, len = {}",
                    forward_strain, real_start, len);
      spdlog::error("real_start + len > seq_view.size(), {} + {} > {}",
                    real_start, len, seq_view.size());
      throw std::out_of_range("end position is out of range of sequence");
    }
    seq_view = seq_view.substr(real_start, len);

    auto qual_view = std::optional<std::string_view>{};
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      qual_view = forward_strain ? this->qual : rev_qual.value();
      assert(qual_view.has_value() && "qual_view is empty");
      qual_view = qual_view->substr(real_start, len);
    }
    return std::make_pair(seq_view, qual_view);
  }


  /**
   * @brief
   *
   * @param start
   * @param end
   * @return auto
   */
  auto add_coverage(const std::size_t start, const std::size_t end) noexcept {
    // spdlog::debug("Add {} - {}", start, end);
    coverage.modify(start, end, 1);
    const double unit_length_ratio = 0.01;

    const int end_threshold = 12;
    const int mid_threshold = 15;

    // if (ready_to_assemble()) {
    //   return;
    // }

    auto can_be_assembled = [&](std::size_t start, std::size_t end) {
      auto head_max_coverage = coverage.get(0, start - 1).max_coverage;
      auto mid_min_coverage = coverage.get(start, end - 1).min_coverage;
      auto tail_max_coverage = coverage.get(end, len() - 1).max_coverage;

      if (head_max_coverage >= end_threshold && 
          tail_max_coverage >= end_threshold &&
          mid_min_coverage >= mid_threshold) {
        // spdlog::debug("H({} - {}) = {}", 0, start - 1, head_max_coverage);
        // spdlog::debug("C({} - {}) = {}", start, end, mid_min_coverage);
        // spdlog::debug("T({} - {}) = {}", end + 1, len() - 1, tail_max_coverage);
        return true;
      }
      return false;
    };

    for (auto i = 1; i <= 2; i++) {
      for (auto j = 1; j <= 2; j++) {
        auto st = len() * (i * unit_length_ratio);
        auto ed = len() * (j * unit_length_ratio);
        ed = len() - ed - 1;
        if (can_be_assembled(st, ed)) {
          state = State::WAIT_FOR_ASSEMBLE;
          return;
        }
      }
    }

    // if coverage is enough, we can set the state to WAIT_FOR_ASSEMBLE
  }

  /**
   * @brief Get the number of fragments for a given read ID.
   *
   * This function returns the number of fragments associated with a specific
   * read ID. If the read ID is not found in the fragment count map, a warning
   * message is logged and 0 is returned.
   *
   * @param read_id The ID of the read.
   * @return The number of fragments for the given read ID.
   */
  auto get_fragment_cnt(std::size_t read_id) noexcept {
    if (!fragment_cnt.contains(read_id)) {
      spdlog::warn("read {} not contains fragments from read {}", this->id,
                   read_id);
      return 0ul;
    }
    return fragment_cnt[read_id];
  }

  /**
   * @brief Get the index range for fragments.
   *
   * This function returns a pair of indices representing the range of fragments
   * to be processed. The function atomically updates the total number of
   * fragments and returns the previous and next values.
   *
   * @param cnt The number of fragments to be processed.
   * @return std::pair<std::size_t, std::size_t> The range of indices.
   */
  auto get_index_range(std::size_t cnt) noexcept {
    auto expected = total_fragments.load();
    auto nxt_value = expected + cnt;
    while (!total_fragments.compare_exchange_weak(expected, nxt_value,
                                                  std::memory_order_acq_rel)) {
      nxt_value = expected + cnt;
    }
    return std::make_pair(expected, nxt_value);
  }

  /**
   * @brief Retrieves the index for a given ID.
   *
   * This function returns the index of `fragments` associated with the given
   * read ID, which is used when adding fragments into read ID. If the ID is not
   * found in the indexes container, a warning message is logged and 0 is
   * returned.
   *
   * @param id The ID for which to retrieve the index.
   * @return The index associated with the given ID, or 0 if the ID is not
   * found.
   */
  auto get_index(std::size_t id) noexcept {
    if (!indexes.contains(id)) {
      spdlog::warn("read {} not contains fragments from read {}", this->id, id);
      return 0ul;
    }
    return indexes[id].fetch_add(1, std::memory_order_acq_rel);
  }

  /**
   * @brief Resizes the fragments vector based on the value of total_fragments.
   *
   * This function resizes the fragments vector to accommodate the specified
   * number of fragments. The size of the vector is determined by the value of
   * the total_fragments atomic variable.
   *
   * @note This function is noexcept.
   */
  auto set_fragments_size() noexcept {
    fragments.resize(total_fragments.load());
  }

  /**
   * @brief Sets the corrected fragment at the specified index.
   *
   * This function sets the corrected fragment at the specified index in the
   * fragments container. The index must be within the valid range of indices.
   *
   * @param idx The index of the fragment to set.
   * @param fragment The corrected fragment to set.
   */
  auto set_corrected_fragment(std::size_t idx, Sequence<std::string> fragment) {
    assert(idx < fragments.size() && "idx out of range");
    if (fragment.left_bound > fragment.right_bound) {
      spdlog::warn("left_bound {} > right_bound {}", fragment.left_bound,
                   fragment.right_bound);
      fragment.left_bound = fragment.right_bound;
    }
    // add_coverage(fragment.left_bound, fragment.right_bound - 1);
    fragments[idx] = std::move(fragment);
  }


  auto save_corrected_fragments(const fs::path& path) {
    auto fout = std::ofstream(path);
    for (auto& f : fragments) {
      if (f.seq.size() >= 0) {
        fout << f.left_bound << ' ' << f.right_bound << ' ' << f.seq << '\n';
      }
    }
  }

  auto add_corrected_fragment(Sequence<std::string> fragment) {
    if (fragment.empty()) {
      return; 
    }
    add_coverage(fragment.left_bound, fragment.right_bound - 1);
    fragments.push_back(std::move(fragment));
  }

  /**
   * @brief Assembles the corrected sequence from the fragments.
   *
   * This function assembles the corrected sequence from the fragments using the
   * ReadAssembler. The assembled sequence is stored in the corrected_seq member
   * variable, and the state is updated to ASSEMBLED.
   *
   * @return None
   */
  auto assemble_corrected_seq() {
    // {
    //   auto path = fs::path(TMP_PATH) / "fragments" / (this->name + ".txt");
    //   auto fout = std::ofstream(path);
    //   for (auto& f : fragments) {
    //     if (f.seq.size() >= 0) {
    //       fout << f.left_bound << ' ' << f.right_bound << ' ' << f.seq <<
    //       '\n';
    //     }
    //   }
    // }
    ReadAssembler assembler(this->seq.size());
    for (auto& f : fragments) {
      assembler.add_seq(std::move(f));
    }
    this->corrected_seq = assembler.assemble();
    if (this->corrected_seq.empty()) {
      state = State::ASSEMBLED_FAILED;
    } else {
      state = State::ASSEMBLED_SUCCESS;
    }
    return !this->corrected_seq.empty();
  }

  auto clear() {
    rc_seq.clear();
    if (rev_qual.has_value()) {
      rev_qual->clear();
    }
    for (auto& w : windows) {
      w.clear();
    }
    windows.clear();
  }

  auto operator<=>(const TestReadWrapper& rhs) const {
    return this->name <=> rhs.name;
  }

  auto check() const {
    auto need_rc = false;
    for (const auto& overlap : overlap_range) {
      if (!overlap.forward_strain) {
        need_rc = true;
      }
    }
    if (need_rc && this->seq.size() != rc_seq.size()) {
      spdlog::error("seq size not equal {} <-> {}", this->seq.size(),
                    rc_seq.size());
      return false;
    }
    if constexpr (std::same_as<R, bio::FastqRecord<R::encoded>>) {
      if (need_rc) {
        assert(rev_qual.has_value() && "rev_qual is empty");
        if (this->qual.size() != rev_qual->size()) {
          spdlog::error("qual size not equal");
          return false;
        }
      }
    }
    return true;
  }

 private:
  // TODO: coverage -> segment tree like data structure
  // ? we may add a data structure here for recording the coverage covered
  // ? by corrected_fragments, if the coverage is enough, we can assemble
  // ? the corrected read without building MSA.

  /* reverse complement sequence */
  std::string rc_seq;

  /* reverse quality sequence */
  std::optional<std::string> rev_qual;

 public:
  /* read id */
  std::size_t id;

  /* flag for need corrected */
  bool need_corrected{true};

  std::atomic<State> state{State::UNCORRECTED};

  /* overlap reads info with this read */
  std::vector<std::size_t> overlap_reads_id;
  std::ranges::subrange<std::vector<Overlap>::iterator> overlap_range;

  /* windows of this read */
  std::vector<Window> windows;

  /* counter of window that finished generating fragments */
  std::atomic_int finished_windows_cnt{0};

  /* counter of fragments of other reads that used to correct this read */
  std::map<std::size_t, std::size_t> fragment_cnt;
  std::map<std::size_t, std::pair<std::size_t, std::size_t>> index_range;
  std::map<std::size_t, std::atomic<std::size_t>> indexes;

  /**
   * corrected sequence fragments from windows of others read, may further
   * assemble to correct version of this read
   */
  std::vector<Sequence<std::string>> fragments;
  std::atomic_uint32_t total_fragments{0};
  segtree coverage{1};

  std::string corrected_seq;
};
