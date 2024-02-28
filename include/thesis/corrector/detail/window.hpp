#pragma once


#include <algorithm>
#include <cassert>
#include <numeric>
#include <optional>
#include <string_view>
#include <vector>
#include <chrono>
#include <random>
#include <ranges>

#include <biovoltron/utility/istring.hpp>
#include <spoa/spoa.hpp>
#include <spdlog/spdlog.h>

#include "thesis/corrector/detail/sequence.hpp"

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

  /**
   * @brief Construct a new Window object
   *
   * @param idx_ index of this window
   * @param start_ start position of this window on raw read sequence
   * @param end_ end position of this window on raw read sequence
   */
  Window(const std::size_t read_id_, const std::size_t idx_,
         const std::size_t start_, const std::size_t end_,
         const std::int64_t seed) {
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
    double PRUNE_LEN_RATIO = 0.95;
  } param;


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
  auto extend(std::size_t extend_len, std::size_t max_len) noexcept {
    start = start > extend_len ? start - extend_len : 0UL;
    end = std::min(end + extend_len, max_len);
  }

  
  /**
   * @brief Returns the depth of the window.
   * 
   * The depth of the window is equal to the number of overlap sequences it contains.
   * 
   * @return The depth of the window.
   */
  auto depth() const noexcept { return overlap_seqs.size(); }


  auto add_sequence(Sequence<>&& seq,
                    const std::pair<std::size_t, std::size_t>& match_pos) {
    if (match_pos.first == match_pos.second || seq.len() == 0) {
      return;
    }
    // TODO: consider about quality filter, low quality -> discard sequence
    overlap_seqs.emplace_back(std::move(seq));
    match_pos_at_query.emplace_back(match_pos);
  }


  auto get_order() noexcept {
    auto order = std::vector<std::size_t>(overlap_seqs.size());
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


  auto build_variation_graph(
    const std::unique_ptr<spoa::AlignmentEngine>& global_aln_engine,
    const std::unique_ptr<spoa::AlignmentEngine>& local_aln_engine,
    const int depth = -1) {
    auto graph = spoa::Graph{};

    auto align_and_push = [&](const auto& seq, const auto& qual, const auto& match_pos) {
      const auto [q_match_st, q_match_ed] = match_pos;
      const auto use_local_aln_ratio = 0.02;
      const auto use_local_aln_len = len() * use_local_aln_ratio;
      auto engine = (spoa::AlignmentEngine*) nullptr;
      if (q_match_st >= use_local_aln_len || q_match_ed <= len() - use_local_aln_len) {
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
    auto order = get_order();
    for (const auto& x : order | std::views::take(depth)) {
      align_and_push(overlap_seqs[x].seq, overlap_seqs[x].qual,
                     match_pos_at_query[x]);
    }
    graph = graph.PruneGraph(len() * param.PRUNE_LEN_RATIO);
    return graph;

    
    // auto order = std::vector<std::size_t>(overlap_seqs.size());
    // std::iota(order.begin(), order.end(), 0);
    // // sort(order.begin(), order.end(), [&](int a, int b) {
    // //   return match_pos_at_query[a].first < match_pos_at_query[b].first;
    // // });
    // std::shuffle(order.begin(), order.end(), std::mt19937{std::random_device{}()});
    // align_and_push(backbone.seq, backbone.qual, std::make_pair(0, len()));
    // for (const auto& x : order) {
    //   align_and_push(overlap_seqs[x].seq, overlap_seqs[x].qual,
    //                  match_pos_at_query[x]);
    // }
    // graph = graph.PruneGraph(len() * param.PRUNE_LEN_RATIO);
    // // graph = graph.PruneGraph(len());

    // return graph;

    // if (read_id < 10) {
    //   auto path =
    //       fs::path(fmt::format("{}/graph2/{}_{}.dot", TMP_PATH, read_id,
    //       idx));
    //   graph.PrintDot(path);
    // }

    // auto path =
    // fs::path(fmt::format("/mnt/ec/ness/yolkee/thesis/tests/graph/{}.dot",
    // idx));
    // graph.PrintDot(path);

    // ? The begin and end may need adjust due to sequencing error
    // ? happen at the begin and end of the read.
    // ? possible solution
    // ? 1. calcuate the coverage for first 10 and last 10 base, select
  }

  auto get_corrected_fragments(
      const std::unique_ptr<spoa::AlignmentEngine>& global_aln_engine,
      const std::unique_ptr<spoa::AlignmentEngine>& local_aln_engine,
      const std::vector<bool> mask,
      const int depth = -1) {
    assert(global_aln_engine != nullptr && local_aln_engine != nullptr);

    if (std::ranges::count(mask, true) == 0) {
      return std::vector<Sequence<std::string>>{};
    }

    auto graph = build_variation_graph(global_aln_engine, local_aln_engine, depth);

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
      //     aln.begin(), std::ranges::find(aln, -1, &std::pair<int, int>::first));
      // // std::ranges::find_last exists in C++23, hell C++ committee
      // const auto tail_not_align = std::ranges::distance(
      //     aln.rbegin(), std::ranges::find(aln.rbegin(), aln.rend(), -1,
      //                                     &std::pair<int, int>::first));
      auto corrected_fragment = graph.DecodeAlignment(aln);
      if (seq.forward_strain == false) {
        corrected_fragment = bio::Codec::rev_comp(corrected_fragment);
      }

      return Sequence<std::string>{
          seq.read_id,                       // read_id
          seq.left_bound,                    // left_bound
          seq.right_bound,                   // right_bound
          std::move(corrected_fragment),     // seq
          std::string{},                     // qual
          seq.forward_strain                 // forward_strain
      };

      // TODO: left and right bound is not corrected
      // (right - left + 1) should equal to ot at least close to
      // corrected_fragment.size()
      // return Sequence<std::string>{
      //     seq.read_id,                       // read_id
      //     seq.left_bound + head_not_align,   // left_bound
      //     seq.right_bound - tail_not_align,  // right_bound
      //     std::move(corrected_fragment),     // seq
      //     std::string{},                     // qual
      //     seq.forward_strain                 // forward_strain
      // };
    };

    auto corrected_fragments = std::vector<Sequence<std::string>>{};
    corrected_fragments.reserve(overlap_seqs.size() + 1);
    corrected_fragments.emplace_back(get_corrected_fragment(backbone));
    for (auto i = 1; i < overlap_seqs.size(); ++i) {
      if (mask[i]) {
        corrected_fragments.emplace_back(get_corrected_fragment(overlap_seqs[i]));
      }
    }
    return corrected_fragments;
  }

  // auto print_graph(const fs::path& path, std::size_t target_read_id) {
  //   auto aln = spoa::Alignment{};
  //   for (auto& seq : overlap_seqs) {
  //     if (seq.read_id != target_read_id) {
  //       continue;
  //     }
  //     aln = spoa::Alignment{};
  //     auto aln_engine = get_global_alignment_engine(500);
  //     aln = aln_engine->Align(seq.seq.data(), seq.seq.size(), graph);
  //     break;
  //   }
  //   graph.PrintDot(path, aln);
  // }

  // auto get_corrected_fragment(
  //     const std::string_view raw_seq) const noexcept {
  //   auto aln = local_aln_engine->Align(raw_seq.data(), raw_seq.size(),
  //   graph); auto corrected_fragment = graph.DecodeAlignment(aln); return
  //   corrected_fragment;
  // };

  // void print_graph_info() const noexcept {
  //   spdlog::debug("read_id = {}", read_id);
  //   spdlog::debug("idx = {}", idx);
  //   spdlog::debug("nodes.size() = {}", graph.nodes().size());
  //   spdlog::debug("edges.size() = {}", graph.edges().size());
  //   spdlog::debug("sequences.size() = {}", graph.sequences().size());
  //   spdlog::debug("rank_to_node.size() = {}", graph.rank_to_node().size());
  //   spdlog::debug("num_codes = {}\n\n", graph.num_codes());
  // }

  void clear() {
    // spdlog::debug("Window {} is destructed", idx);
    match_pos_at_query.clear();
    match_pos_at_query.shrink_to_fit();
    overlap_seqs.clear();
    overlap_seqs.shrink_to_fit();
  }

  /* read id that this window belong */
  std::size_t read_id;

  /* the idx of this window */
  std::size_t idx;

  /* the start and end coordinate of this window on origin read */
  std::size_t start, end;

  /* backbone sequence of this window */
  Sequence<> backbone;

  /* the sequences that overlap with backbone seq */
  std::vector<std::pair<std::size_t, std::size_t>> match_pos_at_query;
  std::vector<Sequence<>> overlap_seqs;

  std::mt19937 random_engine;

  /* debug flag */
  bool debug = false;

  ~Window() { clear(); }
};