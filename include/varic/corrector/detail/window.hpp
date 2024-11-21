#pragma once

#include <algorithm>
#include <boost/asio/detail/posix_thread.hpp>
#include <cassert>
#include <chrono>
#include <numeric>
#include <optional>
#include <random>
#include <ranges>
#include <string_view>
#include <vector>

#include <spdlog/spdlog.h>
#include <biovoltron/utility/istring.hpp>
#include <spoa/spoa.hpp>

#include "varic/corrector/detail/sequence.hpp"

auto get_local_alignment_engine(const std::size_t max_length = 0,
                                const std::size_t match = 5,
                                const std::size_t mismatch = -4,
                                const std::size_t gap = -8,
                                const std::size_t gap_extension = -6) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kSW,  // Smith-Waterman(local alignment)
      match,                     // match (default parameter form SPOA)
      mismatch,                  // mismatch
      gap,                       // gap
      gap_extension              // gap extension
  );
  if (max_length != 0) {
    aln_engine->Prealloc(max_length, 3);
  }
  return aln_engine;
}

auto get_global_alignment_engine(const std::size_t max_length = 0,
                                 const std::size_t match = 5,
                                 const std::size_t mismatch = -4,
                                 const std::size_t gap = -8,
                                 const std::size_t gap_extension = -6) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW,  // Needleman-Wunsch(global alignment)
      match,                     // match (default parameter form SPOA)
      mismatch,                  // mismatch
      gap,                       // gap
      gap_extension              // gap extension
  );
  if (max_length != 0) {
    aln_engine->Prealloc(max_length, 3);
  }
  return aln_engine;
}

auto get_semi_global_alignment_engine(const std::size_t max_length = 0,
                                      const std::size_t match = 5,
                                      const std::size_t mismatch = -4,
                                      const std::size_t gap = -8,
                                      const std::size_t gap_extension = -6) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kOV,  // Needleman-Wunsch(global alignment)
      match,                     // match (default parameter form SPOA)
      mismatch,                  // mismatch
      gap,                       // gap
      gap_extension              // gap extension
  );
  if (max_length != 0) {
    aln_engine->Prealloc(max_length, 3);
  }
  return aln_engine;
}

class Window {
 public:
  Window() = default;

  // delete copy constructor and copy assignment
  Window(const Window&) = delete;
  Window& operator=(const Window&) = delete;

  // default move constructor and move assignment
  Window(Window&&) = default;
  Window& operator=(Window&&) = default;

  using id_t = std::uint32_t;
  using pos_t = std::uint32_t;

  /**
   * @brief Construct a new Window object
   *
   * @param idx_ index of this window
   * @param start_ start position of this window on raw read sequence
   * @param end_ end position of this window on raw read sequence
   */
  Window(const id_t read_id_, const id_t idx_, const pos_t start_,
         const pos_t end_, const std::int64_t seed) {
    read_id = read_id_;
    idx = idx_;
    start = start_;
    end = end_;
    random_engine.seed(seed);
  }

  struct Param {
    const double SHORT_SEQ_RATIO = 0.02;

    /**
     * The length of longest path in pruned graph must greater than this value *
     * len()
     */
  } param;

  /**
   * @brief return the length of this window
   *
   * @return pos_t
   */
  auto len() const noexcept { return end - start; }

  /**
   * @brief extend `start` and `end` by `window_overlap_len` on both side
   * @param max_len the length of raw read sequence
   */
  auto extend(pos_t extend_len, pos_t max_len) noexcept {
    start = start > extend_len ? start - extend_len : 0UL;
    end = std::min(end + extend_len, max_len);
  }

  /**
   * @brief Shrinks the capacity of the internal containers to fit their size.
   *
   * This function reduces the capacity of the `match_pos_at_query` and
   * `overlap_seqs` containers to fit their current size. This can be useful to
   * optimize memory usage when the containers no longer need to hold additional
   * elements.
   *
   * @note This function does not change the size of the containers, only their
   * capacity. If you want to remove excess elements from the containers, use
   * the `clear` function.
   *
   * @note This function is noexcept, meaning it does not throw any exceptions.
   */
  auto shrink_to_fit() noexcept {
    match_pos_at_query.shrink_to_fit();
    overlap_seqs.shrink_to_fit();
  }

  /**
   * @brief Returns the depth of the window.
   *
   * The depth of the window is equal to the number of overlap sequences it
   * contains.
   *
   * @return The depth of the window.
   */
  auto depth() const noexcept { return overlap_seqs.size(); }

  /**
   * Adds a sequence to the window.
   *
   * @param seq The sequence to be added.
   * @param match_pos The position of the match in the sequence.
   *
   * @remarks If the match position range is empty or the sequence length is
   * zero, the sequence will not be added. Consider implementing a quality
   * filter to discard low-quality sequences.
   */
  auto add_sequence(Sequence<>&& seq,
                    const std::pair<pos_t, pos_t>& match_pos) {
    if (match_pos.first == match_pos.second || seq.len() == 0) {
      return;
    }
    // TODO: consider about quality filter, low quality -> discard sequence
    overlap_seqs.emplace_back(std::move(seq));
    match_pos_at_query.emplace_back(match_pos);
  }

  /**
   * @brief Retrieves the order in which to build the multiple sequence
   * alignment (MSA).
   *
   * This function returns a vector of positions representing the order in which
   * to build the MSA. The order is determined based on the length of the
   * sequences in the `overlap_seqs` container.
   *
   * @return A vector of positions representing the order in which to build the
   * MSA.
   */
  auto generate_build_msa_order() noexcept {
    auto order = std::vector<pos_t>(overlap_seqs.size());
    std::iota(order.begin(), order.end(), 0);

    /* seperate full-length sequence and short sequence */
    std::ranges::sort(order, [&](int a, int b) {
      auto a_len = match_pos_at_query[a].second - match_pos_at_query[a].first;
      auto b_len = match_pos_at_query[b].second - match_pos_at_query[b].first;
      return a_len > b_len;
    });

    /* find the seperate point */
    auto sep = std::ranges::find_if(order, [&](int x) {
      auto diff = match_pos_at_query[x].second - match_pos_at_query[x].first;
      return diff < len() * param.SHORT_SEQ_RATIO;
    });

    /* random shuffle two part */
    std::shuffle(order.begin(), sep, random_engine);
    std::shuffle(sep, order.end(), random_engine);
    return order;
  }

  /**
   * Builds a variation graph using alignment engines and other parameters.
   *
   * @param global_aln_engine A unique pointer to the global alignment engine.
   * @param local_aln_engine A unique pointer to the local alignment engine.
   * @param pruned_ratio The ratio used to prune the graph (default: 0.95).
   * @param depth The depth used to limit the number of iterations (default:
   * -1).
   * @return The built variation graph.
   */
  auto build_variation_graph(
      const std::unique_ptr<spoa::AlignmentEngine>& global_aln_engine,
      const std::unique_ptr<spoa::AlignmentEngine>& local_aln_engine,
      const double pruned_ratio, const int depth) {
    auto graph = spoa::Graph{};

    auto align_and_push = [&](const auto& seq, const auto& qual,
                              const auto& match_pos) {
      const auto [q_match_st, q_match_ed] = match_pos;
      const auto use_local_aln_ratio = 0.02;
      const auto use_local_aln_len = len() * use_local_aln_ratio;
      auto engine = (spoa::AlignmentEngine*) nullptr;
      if (q_match_st >= use_local_aln_len ||
          q_match_ed <= len() - use_local_aln_len) {
        engine = local_aln_engine.get();
      } else {
        engine = global_aln_engine.get();
      }
      auto aln = engine->Align(seq.data(), seq.size(), graph);
      if (qual.has_value()) {
        graph.AddAlignment(aln, seq.data(), seq.size(), qual->data(),
                           qual->size());
      } else {
        /* if data has no quality, the weight is set to 1 */
        graph.AddAlignment(aln, seq.data(), seq.size(), 1);
      }
    };

    align_and_push(backbone.seq, backbone.qual, std::make_pair(0, len()));
    /**
     * Actually, the effect of changing order of MSA is not that huge for our
     * method. Use random order is also acceptable.
     */
    auto msa_order = generate_build_msa_order();
    for (const auto& x : msa_order | std::views::take(depth)) {
      align_and_push(overlap_seqs[x].seq, overlap_seqs[x].qual,
                     match_pos_at_query[x]);
    }
    
    // len() -> w
    // pruned_ratio -> p
    graph = graph.PruneGraph(len() * pruned_ratio);
    return graph;
  }

  /**
   * Retrieves the corrected fragments based on the given alignment engines and
   * mask.
   *
   * @param global_aln_engine A unique pointer to the global alignment engine.
   * @param local_aln_engine A unique pointer to the local alignment engine.
   * @param mask A vector of boolean values indicating which fragments need to
   * be corrected and collected.
   * @param pruned_ratio The ratio used to prune the graph (default: 0.95).
   * @param depth The depth of the variation graph (default: -1).
   * @return A vector of corrected fragments.
   */
  auto get_corrected_fragments(
      const std::unique_ptr<spoa::AlignmentEngine>& global_aln_engine,
      const std::unique_ptr<spoa::AlignmentEngine>& local_aln_engine,
      const std::vector<bool> mask, const double pruned_ratio = 0.95,
      const int depth = -1) {
    assert(global_aln_engine != nullptr && local_aln_engine != nullptr);

    auto corrected_fragments_sz = std::ranges::count(mask, true);
    if (corrected_fragments_sz == 0) {
      return std::vector<Sequence<std::string>>{};
    }

    auto graph = build_variation_graph(global_aln_engine, local_aln_engine,
                                       pruned_ratio, depth);

    auto get_corrected_fragment = [&](const Sequence<>& seq) {
      // spoa::Alignment -> std::pair<int, int>
      // first: node index in the graph, -1 mean not align
      auto aln = spoa::Alignment{};
      if (seq.len() < this->len() * 0.98) {
        aln = local_aln_engine->Align(seq.seq.data(), seq.seq.size(), graph);
      } else {
        aln = global_aln_engine->Align(seq.seq.data(), seq.seq.size(), graph);
      }

      // const auto head_not_align = std::ranges::distance(
      //     aln.begin(), std::ranges::find(aln, -1, &std::pair<int,
      //     int>::first));
      // // std::ranges::find_last exists in C++23, hell C++ committee
      // const auto tail_not_align = std::ranges::distance(
      //     aln.rbegin(), std::ranges::find(aln.rbegin(), aln.rend(), -1,
      //                                     &std::pair<int, int>::first));
      auto corrected_fragment = graph.DecodeAlignment(aln);
      if (seq.forward_strain == false) {
        corrected_fragment = bio::Codec::rev_comp(corrected_fragment);
      }

      return Sequence<std::string>{
          seq.read_id,                    // read_id
          seq.left_bound,                 // left_bound
          seq.right_bound,                // right_bound
          std::move(corrected_fragment),  // seq
          std::nullopt,                   // qual
          seq.forward_strain              // forward_strain
      };
    };

    auto corrected_fragments = std::vector<Sequence<std::string>>{};
    corrected_fragments.reserve(corrected_fragments_sz + 1);
    corrected_fragments.emplace_back(get_corrected_fragment(backbone));
    for (auto i = 0; i < overlap_seqs.size(); ++i) {
      if (mask[i]) {
        corrected_fragments.emplace_back(
            get_corrected_fragment(overlap_seqs[i]));
      }
    }
    return corrected_fragments;
  }

  // void print_graph_info() const noexcept {
  //   spdlog::debug("read_id = {}", read_id);
  //   spdlog::debug("idx = {}", idx);
  //   spdlog::debug("nodes.size() = {}", graph.nodes().size());
  //   spdlog::debug("edges.size() = {}", graph.edges().size());
  //   spdlog::debug("sequences.size() = {}", graph.sequences().size());
  //   spdlog::debug("rank_to_node.size() = {}", graph.rank_to_node().size());
  //   spdlog::debug("num_codes = {}\n\n", graph.num_codes());
  // }

  /**
   * @brief Clears the window by removing all stored data.
   *
   * This function clears the window by removing all stored data, including the
   * match positions at the query and the overlap sequences. After calling this
   * function, the window will be empty.
   */
  void clear() {
    // spdlog::debug("Window {} is destructed", idx);
    match_pos_at_query.clear();
    match_pos_at_query.shrink_to_fit();
    overlap_seqs.clear();
    overlap_seqs.shrink_to_fit();
  }

  /* read id that this window belong */
  id_t read_id;

  /* the idx of this window */
  id_t idx;

  /* the start and end coordinate of this window on origin read */
  pos_t start, end;

  /* backbone sequence of this window */
  Sequence<> backbone;

  /* the sequences that overlap with backbone seq */
  std::vector<std::pair<pos_t, pos_t>> match_pos_at_query;
  std::vector<Sequence<>> overlap_seqs;

  std::mt19937 random_engine;

  /* debug flag */
  bool debug = false;

  ~Window() { clear(); }
};