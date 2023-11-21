#pragma once

#include <string_view>

#include "thesis/ds/graph/read_graph.hpp"

class ReadAssembler {

  struct Param {
    const std::size_t MAX_ITERATIONS = 3ul;

    const std::size_t minimum_occurance = 4ul;
    
    std::size_t kmer_size;
    const std::size_t DECREASE_KMER_SIZE = 5ul;
        
    const std::size_t INIT_END_LEN = 100ul;
    const std::size_t INCREASE_END_LEN = 50ul;
    const std::size_t MAX_END_LEN = 500ul;
  } param;

 public:

  ReadAssembler(std::size_t kmer_size,
                std::size_t min_occ,
                std::size_t read_len) : graph(kmer_size, min_occ, read_len) {
    this->param.kmer_size = kmer_size;
  }

  auto add_seq(const std::string_view sequence,
               const std::size_t left_bound,
               const std::size_t right_bound) {
    if (sequence.size() < param.kmer_size) {
      return;
    }
    graph.add_seq(sequence, left_bound);
  }

  auto assemble() {
    // when cannot find proper source and sink
    // 1. decrease kmer size
    // 2. increase end length
    // 3. use backbone sequence as source and sink

    graph.print();

    auto get_sources = [&]() {
      auto end_len = param.INIT_END_LEN;
      while (end_len < param.MAX_END_LEN) {
        auto sources = graph.get_sources(end_len);
        if (sources.size() == 0) {
          spdlog::debug("Head: try {} failed", end_len);
          end_len += param.INCREASE_END_LEN;
          continue;
        }
        return sources;
      }
      return std::vector<ReadGraph::Vertex>{};
    };

    auto get_sinks = [&]() {
      auto end_len = param.INIT_END_LEN;
      while (end_len < param.MAX_END_LEN) {
        auto sinks = graph.get_sinks(end_len);
        if (sinks.size() == 0) {
          spdlog::debug("Tail: try {} failed", end_len);
          end_len += param.INCREASE_END_LEN;
          continue;
        }
        return sinks;
      }
      return std::vector<ReadGraph::Vertex>{};
    };

    auto sources = get_sources();
    auto sinks = get_sinks();

    {
      spdlog::debug("sources.size() = {}", sources.size());
      for (const auto &source : sources) {
        spdlog::debug("source: {}({} - {}) -> {}", graph.g[source].kmer,
                      *graph.g[source].appearances.begin(),
                      *graph.g[source].appearances.rbegin(),
                      graph.g.out_degree(source, false));
      }
      spdlog::debug("sinks size = {}", sinks.size());
      for (const auto& sink : sinks) {
        spdlog::debug("sink: {}({} - {}) -> {}", graph.g[sink].kmer,
                      *graph.g[sink].appearances.begin(),
                      *graph.g[sink].appearances.rbegin(),
                      graph.g.in_degree(sink, false));
      
      }
    }
    if (sources.size() == 0 || sinks.size() == 0) {
      return std::string{};
    }

    for (const auto& source : sources) {
      for (const auto& sink : sinks) {
        auto corrected_read = graph.get_read(source, sink);
        // TODO: add check for correctness
        // TODO: maybe do global or local alignment with raw read
        if (corrected_read.size() != 0) {
          return corrected_read;
        }
      }
    }
    return std::string{};
  }

  auto get_graph_size() const {
    return graph.get_graph_size();
  }

 private:
  std::vector<std::string> seqs;
  ReadGraph graph;
};