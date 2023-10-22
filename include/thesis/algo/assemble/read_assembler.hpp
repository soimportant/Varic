#pragma once

#include <string_view>

#include "thesis/ds/graph/read_graph.hpp"

class ReadAssembler {
 public:
  ReadAssembler(std::size_t kmer_size) : graph(kmer_size, 1) {
  // spdlog::info("ReadAssembler::ReadAssembler()");
}

  auto add_seq(const std::string_view seq, const std::size_t idxL,
               const std::size_t idxR) {
    // spdlog::info("ReadAssembler::add_seq()");
    // TODO:
  }

 private:
  ReadGraph graph;
};