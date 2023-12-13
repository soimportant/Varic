#pragma once

#include "thesis/corrector/detail/sequence.hpp"
#include "thesis/ds/graph/read_graph.hpp"
#include <cmath>
#include <string_view>

class ReadAssembler {

  struct Param {

    /**
     * The parameter that used on finding source and sink
     * at both ends of the read
     */
    /* The length that increase when cannot find source and sink */
    const std::size_t INCREASE_END_LEN = 50ul;
    /* The maximum length for searching, stop searching source and sink when
     * reach this value */
    const std::size_t MAX_END_LEN = 2000ul;

    /**
     * The parameter that used on assembing corrected read
     */

    /* The maximum number of iterations */
    const std::size_t MAX_ITERATIONS = 4ul;
    /* The number of kmer size increase when assemble failed */
    const std::size_t INCREASE_KMER_SIZE = 5ul;

    /**
     * The length ratio that used on checking the result of assembly
     */
    const double ASSEMBLE_FAILED_RATIO = 0.1;
    const double ASSEMBLE_TOO_SHORT_RATIO = 0.9;
    const double ASSEMBLE_MAX_LEN_RATIO = 1.05;
  } param;

  auto identify_source_from_graph(ReadGraph &graph) {
    for (auto start = 0u; start < std::min(read_len, param.MAX_END_LEN);
         start += param.INCREASE_END_LEN) {
      auto end = start + param.INCREASE_END_LEN;
      if (end > read_len) {
        end = read_len;
      }
      auto sources = graph.get_sources(start, end);
      if (sources.size() != 0) {
        return sources;
      }
    }
    return std::vector<ReadGraph::Vertex>{};
  }

  auto identify_sink_from_graph(ReadGraph &graph) {
    for (auto start = 0u; start < std::min(read_len, param.MAX_END_LEN);
         start += param.INCREASE_END_LEN) {
      auto end = start + param.INCREASE_END_LEN;
      if (end > read_len) {
        end = read_len;
      }
      auto sinks = graph.get_sinks(read_len - end, read_len - start);
      if (sinks.size() != 0) {
        return sinks;
      }
    }
    return std::vector<ReadGraph::Vertex>{};
  }

  /**
   * @brief Assembles a read from a sequence of kmers.
   *
   * This function takes a kmer size and a minimum occurrence count, constructs
   * a ReadGraph, and attempts to assemble a read from the graph. It first adds
   * all sequences to the graph, then identifies potential source and sink
   * vertices. If no sources or sinks can be found, it logs a warning and
   * returns an empty string. Otherwise, it attempts to assemble a read from the
   * first source to the first sink. If the assembled read is non-empty, it is
   * returned; otherwise, an empty string is returned.
   *
   * @param kmer_size The size of the kmers to be used in the ReadGraph.
   * @param min_occ The minimum occurrence count for a kmer to be included in
   * the ReadGraph.
   * @return A string representing the assembled read, or an empty string if no
   * read could be assembled.
   */
  auto assemble_read(const std::size_t kmer_size, const std::size_t min_occ) {
    auto graph =
        ReadGraph(kmer_size, min_occ, read_len * param.ASSEMBLE_MAX_LEN_RATIO);
    for (const auto &seq : seqs) {
      graph.add_seq(seq.seq, seq.left_bound);
    }
    auto sources = identify_source_from_graph(graph);
    auto sinks = identify_sink_from_graph(graph);

    // auto sources = std::vector<ReadGraph::Vertex>{};
    // {
    //   for (auto iter = 0u; iter * param.INCREASE_END_LEN < param.MAX_END_LEN;
    //        ++iter) {
    //     auto end_len = iter * param.INCREASE_END_LEN;
    //     if (end_len >= read_len) {
    //       break;
    //     }
    //     sources = graph.get_sources(end_len, end_len +
    //     param.INCREASE_END_LEN); if (sources.size() != 0) {
    //       break;
    //     }
    //   }
    // }
    // auto sinks = std::vector<ReadGraph::Vertex>{};
    // {
    //   for (auto iter = 0u; iter * param.INCREASE_END_LEN <
    //                        std::min(read_len, param.MAX_END_LEN);
    //        ++iter) {
    //     auto end_len = iter * param.INCREASE_END_LEN;
    //     if (end_len + param.INCREASE_END_LEN > read_len) {
    //       break;
    //     }
    //     sinks = graph.get_sinks(read_len - end_len - param.INCREASE_END_LEN,
    //                             read_len - end_len);
    //     if (sinks.size() != 0) {
    //       break;
    //     }
    //   }
    // }

    if (sources.size() == 0 || sinks.size() == 0) {
      spdlog::warn("Cannot find source or sink");
      return std::string{};
    }
    {
      spdlog::debug("Sources size = {}", sources.size());
      for (auto source : sources) {
        spdlog::debug("source = {}({} - {})", graph.g[source].kmer,
                      *graph.g[source].appearances.begin(),
                      *graph.g[source].appearances.rbegin());
      }
      spdlog::debug("Sinks size = {}", sinks.size());
      for (auto sink : sinks) {
        spdlog::debug("sink = {}({} - {})", graph.g[sink].kmer,
                      *graph.g[sink].appearances.begin(),
                      *graph.g[sink].appearances.rbegin());
      }
    }

    auto longest_read = std::string{};
    for (auto source : sources) {
      for (auto sink : sinks) {
        auto read = graph.get_read(source, sink);
        if (read.size() >= read_len * param.ASSEMBLE_TOO_SHORT_RATIO) {
          return read;
        }
        if (read.size() > longest_read.size()) {
          longest_read = read;
        }
      }
    }
    return longest_read;
  }

public:
  ReadAssembler(std::size_t read_len) { this->read_len = read_len; }

  auto add_seq(const std::string_view sequence, const std::size_t left_bound,
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
    auto kmer_size = std::ceil(std::log2(read_len) * 2);
    // TODO: adjust it by coverage of fragments or something else
    auto min_occ = 4;

    for (std::size_t i = 0; i < param.MAX_ITERATIONS; ++i) {
      spdlog::debug("Assemble: kmer_size = {}, min_occ = {}", kmer_size,
                    min_occ);
      auto corrected_read = assemble_read(kmer_size, min_occ);
      if (corrected_read.size() >= read_len * param.ASSEMBLE_TOO_SHORT_RATIO) {
        return corrected_read;
      }
      if (corrected_read.size() > read_len * param.ASSEMBLE_FAILED_RATIO) {
        spdlog::debug("Assemble too short(len = {}), increase kmer size({})",
                      corrected_read.size(), kmer_size);
        kmer_size += param.INCREASE_KMER_SIZE;
      } else {
        spdlog::debug("Assemble failed(len = {}), decrease min_occ({})",
                      corrected_read.size(), min_occ);
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
    //         spdlog::debug("read length too short: {} < {}",
    //         corrected_read.size(), read_len);
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