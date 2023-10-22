#pragma once

#include <numeric>
#include <optional>
#include <spoa/spoa.hpp>
#include <string_view>
#include <vector>

#include "thesis/corrector/detail/sequence.hpp"

// ? maybe semi-global alignment?
auto get_local_alignment_engine(const std::size_t max_length) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kSW,  // Smith-Waterman(local alignment)
      3,                         // match (default parameter form SPOA)
      -5,                        // mismatch
      -4                         // gap
  );
  aln_engine->Prealloc(max_length * 1.5, 5);
  return aln_engine;
};

auto get_global_alignment_engine(const std::size_t max_length) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW,  // Needleman-Wunsch(global alignment)
      3,                         // match (default parameter form SPOA)
      -5,                        // mismatch
      -4                         // gap
  );
  aln_engine->Prealloc(max_length * 1.5, 5);
  return aln_engine;
};

class Window {
 public:
  Window() = default;

  // delete copy constructor and copy assignment
  Window(const Window&) = delete;
  Window& operator=(const Window&) = delete;

  // default move constructor and move assignment
  Window(Window&&) = default;
  Window& operator=(Window&&) = default;

  /**
   * @brief Construct a new Window object
   *
   * @param idx_ index of this window
   * @param start_ start position of this window on raw read sequence
   * @param end_ end position of this window on raw read sequence
   */
  Window(const std::size_t idx_, const std::size_t start_,
         const std::size_t end_) {
    idx = idx_;
    start = start_;
    end = end_;
  }

  /**
   * @brief return the length of this window
   *
   * @return std::size_t
   */
  auto len() const noexcept { return end - start; }

  /**
   * @brief extend `start` and `end` by `window_overlap_len` on both side
   * @param max_len the length of raw read sequence
   */
  auto extend(const std::size_t max_len) noexcept {
    start = start > window_overlap_len ? start - window_overlap_len : 0UL;
    end = std::min(end + window_overlap_len, max_len);
  }

  // basic length of a window
  static const auto window_len = 500ul;
  // overlap length of adjanency window. Therefore, typeical length of a
  // single window should be window_len + 2 * window_overlap_len
  static const auto window_overlap_len = 25ul;

  // minimum length of overlap
  static const auto overlap_min_len = window_len - window_overlap_len;

  // auto print() const noexcept {
  //   spdlog::debug("idxL = {}, idxR = {}", idxL, idxR);
  //   spdlog::debug("backbone_seq.size() = {}", backbone_seq.size());
  //   auto sum = std::accumulate(overlap_seqs.begin(), overlap_seqs.end(),
  //   0ul,
  //                              [](auto a, auto& b) { return a + b.size();
  //                              });
  //   spdlog::debug("average overlap_seqs len = {}", sum /
  //   overlap_seqs.size()); spdlog::debug("get prune len = {}",
  //   get_prune_len());
  // }

  auto build_variation_graph(
      std::shared_ptr<spoa::AlignmentEngine> aln_engine) {
    assert(aln_engine != nullptr);
    auto align_and_push = [&](const std::string_view seq,
                              const std::optional<std::string_view>& qual) {
      auto aln = aln_engine->Align(seq.data(), seq.size(), graph);
      if (qual.has_value()) {
        graph.AddAlignment(aln, seq.data(), seq.size(), qual.value().data(),
                           qual.value().size());
      } else {
        graph.AddAlignment(aln, seq.data(), seq.size(), 1);
      }
    };

    // TODO: random_view

    auto order = std::vector<std::size_t>(overlap_seqs.size());
    std::iota(order.begin(), order.end(), 0);
    // std::random_shuffle(order.begin(), order.end());
    align_and_push(backbone.seq, backbone.qual);
    for (const auto& now : order) {
      align_and_push(overlap_seqs[now].seq, overlap_seqs[now].qual);
    }
    // ? The begin and end may need adjust due to sequencing error
    // ? happen at the begin and end of the read.
    // ? possible solution
    // ? 1. calcuate the coverage for first 10 and last 10 base, select
  }

  /* the idx of this window */
  std::size_t idx;

  /* the range covered by this window */
  std::size_t start, end;

  /* backbone sequence of this window */
  Sequence backbone;

  /* the sequence that overlap with backbone seq */
  std::vector<Sequence> overlap_seqs;

 private:
  /* the variation graph */
  spoa::Graph graph;
};