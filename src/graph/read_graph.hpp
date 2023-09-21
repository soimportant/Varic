#pragma once

#include <map>
#include <set>
#include <string>

#include <biovoltron/algo/assemble/graph/graph_wrapper.hpp>
#include <spdlog/spdlog.h>

struct ReadGraph {
private:
  struct VertexProperty {
    std::string kmer;
    std::size_t idxL;
    std::size_t idxR;
  };

  struct EdgeProperty {
    std::size_t count = 0;
  };

  using Graph = bio::GraphWrapper<VertexProperty, EdgeProperty>;
  using Vertex = Graph::Vertex;
  using Edge = Graph::Edge;
  using Path = std::vector<Vertex>;

  Graph g;
  // ? the usage of dup_kmers
  std::set<std::string> dup_kmers;
  // ? the usage of unique_kmers
  std::map<std::string, Vertex> unique_kmers;

  std::size_t kmer_size;
  std::size_t minimum_occurance;

  auto build_dup_kmers(const std::string_view seq) {
    std::set<std::string_view> kmers;
    for (auto i = 0u; i + kmer_size <= seq.size(); i++) {
      auto subseq = seq.substr(i, kmer_size);
      if (kmers.contains(subseq)) {
        dup_kmers.insert(std::string(subseq));
      }
      kmers.insert(subseq);
    }
  }

  auto create_vertex(const std::string kmer) {

    // auto [iter, success] = unique_kmers.insert_or_assign(kmer, Vertex{});

    if (unique_kmers.contains(kmer)) {
      return unique_kmers[kmer];
    }
    const auto v = g.create_vertex();
    g[v].kmer = kmer;
    unique_kmers[kmer] = v;
    return v;
  }

  auto create_edge(Vertex u, Vertex v) {
    const auto e = g.create_edge(u, v);
    g[e].count += 1;
    return e;
  }

  auto extend_chain(const Vertex& u, const std::string kmer) {
    for (const auto e : g.out_edges(u, false)) {
      const auto v = g.target(e);
      if (g[v].kmer == kmer) {
        g[e].count += 1;
        return v;
      }
    }
    const auto v = create_vertex(kmer);
    create_edge(u, v);
    return v;
  }

  auto get_vertex(const std::string& kmer) {
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) {
      return it->second;
    }
    return create_vertex(kmer);
  }

  auto concat_vertices(const Path& path) {
    auto u = path[0];
    auto seq = g[u].kmer;
    for (auto i = 1u; i < path.size(); i++) {
      auto v = path[i];
      seq += g[v].kmer.back();
    }
    return seq;
  }

public:
  ReadGraph(const std::size_t kmer_size, const std::size_t minimum_occurance)
      : kmer_size(kmer_size), minimum_occurance(minimum_occurance) {
    g.set_edge_filter([this](const Edge& e) {
      return g[e].count >= this->minimum_occurance;
    });
  }

  /**
   * @brief add one sequence to the graph
   * @param sequence
   * @note provide API for adding only one sequence for reducing memory
   * consumption
   */
  void add_seq(const std::string& sequence) {
    // auto sequence = std::string_view(sequence);
    auto v = get_vertex(sequence.substr(0, kmer_size));
    for (auto i = 1u; i + kmer_size <= sequence.size(); i++) {
      v = extend_chain(v, sequence.substr(i, kmer_size));
    }
  }

  /**
   * @brief add one sequence to the graph
   * @param sequence
   * @note provide API for adding only one sequence for reducing memory
   * consumption
   */
  void add_seq(const std::string& sequence, std::size_t idxL, std::size_t idxR) {
    // auto sequence = std::string_view(sequence);
    auto v = get_vertex(sequence.substr(0, kmer_size));
    for (auto i = 1u; i + kmer_size <= sequence.size(); i++) {
      v = extend_chain(v, sequence.substr(i, kmer_size));
    }
  }

  auto build(const std::vector<std::string>& sequences) {
    // for (auto& seq : sequences) {
    //   build_dup_kmers(seq);
    // }
    for (auto& seq : sequences) {
      if (seq.size() >= kmer_size) {
        add_seq(seq);
      }
    }
  }

  void print() {
    auto sources = g.get_sources();
    auto sinks = g.get_sinks();
    auto nodes_size = std::ranges::ssize(g.vertices());

    // for (auto v : g.vertices() | std::views::take(2)) {
    //   spdlog::info("vertex: {}", g[v].kmer);
    //   for (auto e : g.out_edges(v, false)) {
    //     auto u = g.target(e);
    //     spdlog::info("{} -> {}: {}", g[v].kmer, g[u].kmer, g[e].count);
    //   }
    // }

    spdlog::info("sources: {}", sources.size());
    spdlog::info("sinks: {}", sinks.size());
    /* select front kmer as sources */
    /* select the last kmer as sinks */
    // for (auto &source : sources) {
    //   spdlog::info("source: {}", g[source].kmer);
    // }
    // for (auto &sink : sinks) {
    //   spdlog::info("sink: {}", g[sink].kmer);
    // }

  }
};