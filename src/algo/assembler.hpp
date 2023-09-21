#pragma once

#include <vector>
#include <string>

#include <spdlog/spdlog.h>
#include "../graph/read_graph.hpp"


class ReadAssembler {
public:
  ReadAssembler(std::size_t kmer_size) : graph(kmer_size, 1) {
    
  }

  auto add_seq(const std::string_view seq,
               const std::size_t idxL,
               const std::size_t idxR) {
    // graph.add_seq(seq, idxL, idxR);
  }

private:
  ReadGraph graph;
};