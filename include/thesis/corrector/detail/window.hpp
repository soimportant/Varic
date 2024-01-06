#pragma once

#include <algorithm>
#include <cassert>
#include <numeric>
#include <optional>
#include <string_view>
#include <vector>

#include <spdlog/spdlog.h>
#include <spoa/spoa.hpp>

#include "thesis/corrector/detail/sequence.hpp"

// ? maybe semi-global alignment?
auto get_local_alignment_engine(const std::size_t max_length = 0) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kSW, // Smith-Waterman(local alignment)
      5,                        // match (default parameter form SPOA)
      -4,                       // mismatch
      -8,                       // gap
      -6                        // gap extension
  );
  // 5 -4 -8 -6
  if (max_length != 0) {
    aln_engine->Prealloc(max_length, 5);
  }
  return aln_engine;
};

auto get_global_alignment_engine(const std::size_t max_length = 0) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, // Needleman-Wunsch(global alignment)
      5,                        // match (default parameter form SPOA)
      -4,                       // mismatch
      -8,                       // gap
      -6                        // gap extension
  );
  if (max_length != 0) {
    aln_engine->Prealloc(max_length, 5);
  }
  return aln_engine;
};

class Window {
public:
  Window() = default;

  // delete copy constructor and copy assignment
  Window(const Window &) = delete;
  Window &operator=(const Window &) = delete;

  // default move constructor and move assignment
  Window(Window &&) = default;
  Window &operator=(Window &&) = default;

  /**
   * @brief Construct a new Window object
   *
   * @param idx_ index of this window
   * @param start_ start position of this window on raw read sequence
   * @param end_ end position of this window on raw read sequence
   */
  Window(const std::size_t read_id_, const std::size_t idx_,
         const std::size_t start_, const std::size_t end_) {
    read_id = read_id_;
    idx = idx_;
    start = start_;
    end = end_;
  }

  struct Param {
    /**
     * The length of longest path in pruned graph must greater than this value *
     * len()
     */
    const double PRUNE_LEN_RATIO = 0.98;
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

  auto add_sequence(Sequence<> seq,
                    const std::pair<std::size_t, std::size_t> &match_pos) {
    if (match_pos.first == match_pos.second || seq.len() == 0) {
      return;
    }
    // TODO: consider about quality filter, low quality -> discard sequence
    overlap_seqs.emplace_back(std::move(seq));
    match_pos_at_query.emplace_back(match_pos);
  }

  auto
  build_variation_graph(std::shared_ptr<spoa::AlignmentEngine> aln_engine) {
    assert(aln_engine != nullptr);
    auto align_and_push =
        [&](const std::string_view seq,
            const std::optional<std::string_view> &qual,
            const std::pair<std::size_t, std::size_t> query_match_range) {
          const auto [q_match_st, q_match_ed] = query_match_range;
          auto aln = spoa::Alignment{};
          const auto align_on_subgraph_len = len() * 0.02;
          if (query_match_range.first >= align_on_subgraph_len &&
              query_match_range.second <= len() - align_on_subgraph_len) {
            auto node_mapping = std::vector<const spoa::Graph::Node *>{};
            // spdlog::debug("q_match_st = {}, q_match_ed = {}", q_match_st,
            // q_match_ed);
            auto subgraph =
                graph.Subgraph(q_match_st, q_match_ed - 1, &node_mapping);
            aln = aln_engine->Align(seq.data(), seq.size(), subgraph);
            subgraph.UpdateAlignment(node_mapping, &aln);
          } else {
            aln = aln_engine->Align(seq.data(), seq.size(), graph);
          }
          if (qual.has_value()) {
            graph.AddAlignment(aln, seq.data(), seq.size(), qual->data(),
                               qual->size());
          } else {
            /* if data has no quality, the weight is set to 1 */
            graph.AddAlignment(aln, seq.data(), seq.size(), 1);
          }
        };

    /* TODO: the order is used by vechat, but can try it randomly */
    auto order = std::vector<std::size_t>(overlap_seqs.size());
    std::iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
      return match_pos_at_query[a].first < match_pos_at_query[b].first;
    });
    align_and_push(backbone.seq, backbone.qual, {0, len()});
    for (const auto &x : order) {
      align_and_push(overlap_seqs[x].seq, overlap_seqs[x].qual,
                     match_pos_at_query[x]);
    }
    // graph = graph.PruneGraph(len() * param.PRUNE_LEN_RATIO);
    graph = graph.PruneGraph(len());

    // if (read_id < 10) {
    //   auto path =
    //       fs::path(fmt::format("{}/graph2/{}_{}.dot", TMP_PATH, read_id, idx));
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
      std::shared_ptr<spoa::AlignmentEngine> global_aln_engine) {
    assert(global_aln_engine != nullptr);
    auto local_aln_engine = get_local_alignment_engine();

    auto get_corrected_sequence = [&](const Sequence<> &seq) {
      // spoa::Alignment -> std::pair<int, int>
      // first: node index in the graph, -1 mean not align
      auto aln = spoa::Alignment{};
      if (seq.len() < this->len() * 0.8) {
        aln = local_aln_engine->Align(seq.seq.data(), seq.seq.size(), graph);
      } else {
        aln = global_aln_engine->Align(seq.seq.data(), seq.seq.size(), graph);
      }

      const auto head_not_align = std::ranges::distance(
        aln.begin(),
        std::ranges::find(aln, -1, &std::pair<int, int>::first)
      );
      // std::ranges::find_last exists in C++23, hell C++ committee
      const auto tail_not_align = std::ranges::distance(
        aln.rbegin(),
        std::ranges::find(aln.rbegin(), aln.rend(), -1, 
                          &std::pair<int, int>::first)
      );
      auto corrected_fragment = graph.DecodeAlignment(aln);
      return Sequence<std::string>{
          seq.read_id,                      // read_id
          seq.left_bound + head_not_align,  // left_bound
          seq.right_bound - tail_not_align, // right_bound
          std::move(corrected_fragment),    // seq
          std::string{},                    // qual
          seq.forward_strain                // forward_strain
      };
    };

    auto order = std::vector<std::size_t>(overlap_seqs.size());
    std::iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
      return match_pos_at_query[a].first < match_pos_at_query[b].first;
    });

    auto corrected_sequences = std::vector<Sequence<std::string>>{};
    corrected_sequences.reserve(overlap_seqs.size() + 1);
    corrected_sequences.emplace_back(get_corrected_sequence(backbone));
    for (const auto &x : order) {
      corrected_sequences.emplace_back(get_corrected_sequence(overlap_seqs[x]));
    }
    return corrected_sequences;
  }

  // auto get_corrected_fragment(
  //     const std::string_view raw_seq) const noexcept {
  //   auto aln = local_aln_engine->Align(raw_seq.data(), raw_seq.size(),
  //   graph); auto corrected_fragment = graph.DecodeAlignment(aln); return
  //   corrected_fragment;
  // };

  void print_graph_info() const noexcept {
    spdlog::debug("read_id = {}", read_id);
    spdlog::debug("idx = {}", idx);
    spdlog::debug("nodes.size() = {}", graph.nodes().size());
    spdlog::debug("edges.size() = {}", graph.edges().size());
    spdlog::debug("sequences.size() = {}", graph.sequences().size());
    spdlog::debug("rank_to_node.size() = {}", graph.rank_to_node().size());
    spdlog::debug("num_codes = {}\n\n", graph.num_codes());
  }

  void clear() {
    // spdlog::debug("Window {} is destructed", idx);
    match_pos_at_query.clear();
    overlap_seqs.clear();
    graph.Clear();
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

  ~Window() { clear(); }

private:
  /* the variation graph */
  spoa::Graph graph;
};