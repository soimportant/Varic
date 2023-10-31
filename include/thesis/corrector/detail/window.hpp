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
      // -6                           // gap extension
  );
  // 5 -4 -8 -6
  aln_engine->Prealloc(max_length * 5, 5);
  return aln_engine;
};

auto get_global_alignment_engine(const std::size_t max_length) {
  auto aln_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW,  // Needleman-Wunsch(global alignment)
      5,                         // match (default parameter form SPOA)
      -4,                        // mismatch
      -8                         // gap
  );
  aln_engine->Prealloc(max_length * 5, 5);
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
  Window(const std::size_t read_id_, const std::size_t idx_,
         const std::size_t start_, const std::size_t end_) {
    read_id = read_id_;
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
  auto extend(std::size_t extend_len, std::size_t max_len) noexcept {
    start = start > extend_len ? start - extend_len : 0UL;
    end = std::min(end + extend_len, max_len);
  }

  auto add_sequence(const Sequence& seq,
                    const std::pair<std::size_t, std::size_t>& match_pos) {
    if (match_pos.first == match_pos.second) {
      return;
    }
    if (seq.len() == 0) {
      return;
    }

    overlap_seqs.emplace_back(std::move(seq));
    match_pos_at_query.emplace_back(match_pos);
  }

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
    auto align_and_push =
        [&](const std::string_view seq,
            const std::optional<std::string_view>& qual,
            const std::pair<std::size_t, std::size_t> query_match_range) {
      const auto& [q_match_st, q_match_ed] = query_match_range;
      auto aln = spoa::Alignment{};
      const auto align_on_subgraph_len = len() * 0.02;
      if (query_match_range.first >= align_on_subgraph_len &&
          query_match_range.second <= len() - align_on_subgraph_len) {
        auto node_mapping = std::vector<const spoa::Graph::Node*>{};
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
        graph.AddAlignment(aln, seq.data(), seq.size(), qual.value().data(),
                           qual.value().size());
      } else {
        /* if data has no quality, the weight is set to 1 */
        graph.AddAlignment(aln, seq.data(), seq.size(), 1);
      }
    };

    auto order = std::vector<std::size_t>(overlap_seqs.size());
    std::iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
      return match_pos_at_query[a].first < match_pos_at_query[b].first;
    });
    align_and_push(backbone.seq, backbone.qual, {0, len()});
    for (const auto& x : order) {
      align_and_push(overlap_seqs[x].seq, overlap_seqs[x].qual,
                     match_pos_at_query[x]);
    }

    // auto msa = graph.GenerateMultipleSequenceAlignment();
    // auto out_path =
    // fs::path(fmt::format("/mnt/ec/ness/yolkee/thesis/tests/tmp/msa/window_{}.txt",
    // idx)); std::ofstream out(out_path); out << msa[0].size() << std::endl;
    // for (auto &s : msa) {
    //   out << s << std::endl;
    // }

    // ? The begin and end may need adjust due to sequencing error
    // ? happen at the begin and end of the read.
    // ? possible solution
    // ? 1. calcuate the coverage for first 10 and last 10 base, select

    graph.Clear();
  }

  void print_graph_info() {
    spdlog::debug("read_id = {}", read_id);
    spdlog::debug("idx = {}", idx);
    spdlog::debug("nodes.size() = {}", graph.nodes().size());
    spdlog::debug("edges.size() = {}", graph.edges().size());
    spdlog::debug("sequences.size() = {}", graph.sequences().size());
    spdlog::debug("rank_to_node.size() = {}", graph.rank_to_node().size());
    spdlog::debug("num_codes = {}", graph.num_codes());
  }

  /* read id that this window belong */
  std::size_t read_id;

  /* the idx of this window */
  std::size_t idx;

  /* the start and end coordinate of this window on origin read */
  std::size_t start, end;

  /* backbone sequence of this window */
  Sequence backbone;

  /* the sequence that overlap with backbone seq */
  std::vector<std::pair<std::size_t, std::size_t>> match_pos_at_query;
  std::vector<Sequence> overlap_seqs;

 private:
  /* the variation graph */
  spoa::Graph graph;
};