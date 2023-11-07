#pragma once

#include <string_view>

#include "thesis/ds/graph/read_graph.hpp"

class ReadAssembler {
 public:
  ReadAssembler() = default;

  ReadAssembler(std::size_t kmer_size) : graph(kmer_size, 1) {
    this->kmer_size = kmer_size;
  }

  auto add_seq(const std::string& seq, const std::size_t idxL,
               const std::size_t idxR) {
    // TODO: if assembley result is not good, add position information for
    // TODO: making it more precise
    if (seq.size() <= kmer_size) {
      return;
    }
    graph.add_seq(seq, idxL, idxR);
  }

  auto get_graph_size() const {
    return graph.get_graph_size();
  }

 private:
  std::size_t kmer_size;
  ReadGraph graph;
};