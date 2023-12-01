#pragma once

#include <string_view>

#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/ds/graph/read_graph.hpp"


class ReadAssembler {

  struct Param {
    const std::size_t MAX_ITERATIONS = 3ul;
    const std::size_t MIN_OCC = 4ul;    
    const std::size_t INCREASE_KMER_SIZE = 5ul;
        
    const std::size_t INIT_END_LEN = 100ul;
    const std::size_t INCREASE_END_LEN = 50ul;
    const std::size_t MAX_END_LEN = 500ul;

    const double ASSEMBLE_FAILED_RATIO = 0.1;
    const double ASSEMBLE_TOO_SHORT_RATIO = 0.9;
    const double ASSEMBLE_MAX_LEN_RATIO = 1.05;
  } param;


  auto assemble_read(std::size_t kmer_size,
                std::size_t min_occ) {
    auto graph = ReadGraph(kmer_size, min_occ);
    for (const auto& seq : seqs) {
      graph.add_seq(seq.seq, seq.left_bound);
    }
    auto sources = std::vector<ReadGraph::Vertex>{};
    {
      auto end_len = param.INIT_END_LEN;
      while (end_len < param.MAX_END_LEN) {
        sources = graph.get_sources(end_len);
        if (sources.size() != 0) {
          break;
        }
        spdlog::debug("Head: try {} failed", end_len);
        end_len += param.INCREASE_END_LEN;
      }
    }
    auto sinks = std::vector<ReadGraph::Vertex>{};
    {
      auto end_len = param.INIT_END_LEN;
      while (end_len < param.MAX_END_LEN) {
        sinks = graph.get_sinks(read_len - end_len);
        if (sinks.size() != 0) {
          break;
        }
        spdlog::debug("Tail: try {} failed", end_len);
        end_len += param.INCREASE_END_LEN;
      }
    }
    
    if (sources.size() == 0 || sinks.size() == 0) {
      return std::string{};
    }
    {
      spdlog::debug("Sources size = {}", sources.size());
      for (auto source : sources) {
        spdlog::debug("source = {}({} - {})", graph.g[source].kmer, *graph.g[source].appearances.begin(), *graph.g[source].appearances.rbegin());
      }
      spdlog::debug("Sinks size = {}", sinks.size());
      for (auto sink : sinks) {
        spdlog::debug("sink = {}({} - {})", graph.g[sink].kmer, *graph.g[sink].appearances.begin(), *graph.g[sink].appearances.rbegin());
      }
    }
    // for (const auto& source : sources) {
    //   for (const auto& sink : sinks) {
    //     auto corrected_read = graph.get_read(source, sink, read_len * param.ASSEMBLE_MAX_LEN_RATIO);
    //     spdlog::debug("corrected_read = {}", corrected_read.size());
    //   }
    // }

    auto source = sources[0];
    auto sink = sinks[0];
    auto corrected_read = graph.get_read(source, sink, read_len * param.ASSEMBLE_MAX_LEN_RATIO);
    if (corrected_read.size() != 0) {
      return corrected_read;
    }
    return std::string{};
  }

 public:

  ReadAssembler(std::size_t read_len) {
    this->read_len = read_len;
  }

  auto add_seq(const std::string_view sequence,
               const std::size_t left_bound,
               const std::size_t right_bound) {
    seqs.emplace_back(Sequence<std::string>{
      .read_id = -1,
      .left_bound = left_bound,
      .right_bound = right_bound,
      .seq = std::string(sequence),
      .qual = std::nullopt,
      .forward_strain = true,
    });
  }

  auto assemble() {
    
    // TODO: can be higher
    auto kmer_size = std::ceil(std::log2(read_len) * 3);
    auto min_occ = param.MIN_OCC;

    for (std::size_t i = 0; i < param.MAX_ITERATIONS; ++i) {
      spdlog::debug("Assemble: kmer_size = {}, min_occ = {}", kmer_size, min_occ);
      auto corrected_read = assemble_read(kmer_size, min_occ);
      if (corrected_read.size() >= read_len * param.ASSEMBLE_TOO_SHORT_RATIO) {
        return corrected_read;
      }
      if (corrected_read.size() > read_len * param.ASSEMBLE_FAILED_RATIO) {
        spdlog::debug("Assemble too short(len = {}), increase kmer size({})", corrected_read.size(), kmer_size);
        kmer_size += param.INCREASE_KMER_SIZE;
        min_occ -= 1;
      } else {
        spdlog::debug("Assemble failed(len = {}), decrease min_occ({})", corrected_read.size(), min_occ);
        kmer_size += param.INCREASE_KMER_SIZE;
      }
    }
    return std::string{};

    // when cannot find proper source and sink
    // 1. decrease kmer size
    // 2. increase end length
    // 3. use backbone sequence as source and sink

    // graph.print();

    // auto get_sources = [&]() {
    //   auto end_len = param.INIT_END_LEN;
    //   while (end_len < param.MAX_END_LEN) {
    //     auto sources = graph.get_sources(end_len);
    //     if (sources.size() == 0) {
    //       spdlog::debug("Head: try {} failed", end_len);
    //       end_len += param.INCREASE_END_LEN;
    //       continue;
    //     }
    //     return sources;
    //   }
    //   return std::vector<ReadGraph::Vertex>{};
    // };

    // auto get_sinks = [&]() {
    //   auto end_len = param.INIT_END_LEN;
    //   while (end_len < param.MAX_END_LEN) {
    //     auto sinks = graph.get_sinks(end_len);
    //     if (sinks.size() == 0) {
    //       spdlog::debug("Tail: try {} failed", end_len);
    //       end_len += param.INCREASE_END_LEN;
    //       continue;
    //     }
    //     return sinks;
    //   }
    //   return std::vector<ReadGraph::Vertex>{};
    // };

    // auto sources = get_sources();
    // auto sinks = get_sinks();
    // if (sources.size() == 0 || sinks.size() == 0) {
    //   return std::string{};
    // }

    // for (const auto& source : sources) {
    //   for (const auto& sink : sinks) {
    //     auto corrected_read = graph.get_read(source, sink);
    //     if (corrected_read.size() != 0) {
    //       if (corrected_read.size() < read_len * 0.93) {
    //         spdlog::debug("read length too short: {} < {}", corrected_read.size(), read_len);
    //       }
    //       return corrected_read;
    //     }
    //   }
    // }
    // return std::string{};
  }

 private:
  std::size_t read_len;
  std::vector<Sequence<std::string>> seqs;
};