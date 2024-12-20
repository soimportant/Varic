#pragma once

#include <cmath>
#include <string_view>

#include "varic/corrector/detail/sequence.hpp"
#include "varic/ds/graph/read_graph.hpp"

class ReadAssembler {
  enum State {
    // INIT,
    ASSEMBLE_FAILED,
    ASSEMBLE_TOO_SHORT,
    ASSEMBLE_SUCCESS,
    ASSEMBLE_TOO_LONG,
  };

  struct Param {
    /**
     * Parameters that used on identifying source and sink
     */
    // The ratio of searching length related to `read_len` on both end
    const double INCREASE_END_LEN_RATIO = 0.02;

    // The maximum searching length ratio
    const double MAX_END_LEN_RATIO = 0.12;

    // Maximum occurances at searching source and sink, *User parameter*
    const std::size_t MAX_END_OCC = 6ul;
    const std::size_t MIN_END_OCC = 3ul;

    /**
     * The parameter that used on assembing corrected read
     */
    // The maximum number of iterations
    const std::size_t MAX_ITERATIONS = 4ul;
    // The number of kmer size increase when assemble failed
    const std::size_t INCREASE_KMER_SIZE = 5ul;

    /**
     * The length ratio that used on checking the result of assembly
     */
    const double ASSEMBLE_FAILED_RATIO = 0.1;
    const double ASSEMBLE_TOO_SHORT_RATIO = 0.9;
    const double ASSEMBLE_MAX_LEN_RATIO = 1;
  } param;

  /**
   * Identifies the source vertices in the given ReadGraph based on the
   * specified criteria.
   *
   * @param graph The ReadGraph to analyze.
   * @return A vector of source vertices satisfying the criteria, or an empty
   * vector if no sink vertices are found.
   */
  auto identify_source_from_graph(ReadGraph& graph) {
    for (auto occ = param.MAX_END_OCC; occ >= param.MIN_END_OCC; --occ) {
      for (auto ratio = 0.0; ratio < param.MAX_END_LEN_RATIO;
           ratio += param.INCREASE_END_LEN_RATIO) {
        auto start = std::floor(read_len * ratio);
        auto end = std::ceil(read_len * (ratio + param.INCREASE_END_LEN_RATIO));
        /* get sources from (start, end) with occurance = `occ` */
        auto sources = graph.get_sources(start, end, occ);
        if (sources.size() != 0) {
          return sources;
        }
      }
    }
    return std::vector<ReadGraph::Vertex>{};
  }

  /**
   * Identifies the sink vertices from the given ReadGraph based on the
   * specified criteria.
   *
   * @param graph The ReadGraph to identify the sink vertices from.
   * @return A vector of sink vertices satisfying the criteria, or an empty
   * vector if no sink vertices are found.
   */
  auto identify_sink_from_graph(ReadGraph& graph) {
    for (auto weight = param.MAX_END_OCC; weight >= param.MIN_END_OCC;
         --weight) {
      for (auto ratio = 0.0; ratio < param.MAX_END_LEN_RATIO;
           ratio += param.INCREASE_END_LEN_RATIO) {
        auto start =
            std::ceil(read_len * (1 - (ratio + param.INCREASE_END_LEN_RATIO)));
        auto end = std::floor(read_len * (1 - ratio));
        auto sinks = graph.get_sinks(start, end, weight);
        if (sinks.size() != 0) {
          return sinks;
        }
      }
    }
    return std::vector<ReadGraph::Vertex>{};
  }

  /**
   * Determines the state of a read based on its length.
   *
   * @param read The input read to be evaluated.
   * @return The state of the read (ASSEMBLE_TOO_LONG, ASSEMBLE_SUCCESS,
   * ASSEMBLE_TOO_SHORT, or ASSEMBLE_FAILED).
   */
  auto get_read_state(std::string_view read) {
    if (read.size() >= read_len * param.ASSEMBLE_MAX_LEN_RATIO) {
      return State::ASSEMBLE_TOO_LONG;
    }
    if (read.size() >= read_len * param.ASSEMBLE_TOO_SHORT_RATIO) {
      return State::ASSEMBLE_SUCCESS;
    }
    if (read.size() > read_len * param.ASSEMBLE_FAILED_RATIO) {
      return State::ASSEMBLE_TOO_SHORT;
    }
    return State::ASSEMBLE_FAILED;
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
  auto assemble_read(const std::size_t kmer_size) {
    auto graph = ReadGraph(kmer_size, read_len * param.ASSEMBLE_MAX_LEN_RATIO);
    for (const auto& seq : sequences) {
      graph.add_seq(seq.seq, seq.left_bound);
    }
    auto sources = identify_source_from_graph(graph);
    auto sinks = identify_sink_from_graph(graph);

    if (sources.size() == 0 || sinks.size() == 0) {
      // spdlog::warn("Cannot find source or sink");
      return std::string{};
    }
    // {
    //   spdlog::debug("Sources size = {}", sources.size());
    //   for (auto source : sources) {
    //     spdlog::debug("source = {}({} - {})", graph.g[source].kmer,
    //                   *graph.g[source].appearances.begin(),
    //                   *graph.g[source].appearances.rbegin());
    //   }
    //   spdlog::debug("Sinks size = {}", sinks.size());
    //   for (auto sink : sinks) {
    //     spdlog::debug("sink = {}({} - {})", graph.g[sink].kmer,
    //                   *graph.g[sink].appearances.begin(),
    //                   *graph.g[sink].appearances.rbegin());
    //   }
    // }

    auto longest_read = std::string{};
    for (auto source : sources | std::views::take(2)) {
      for (auto sink : sinks | std::views::take(2)) {
        auto forward_read = graph.find_read<false>(source, sink);
        switch (get_read_state(forward_read)) {
          case State::ASSEMBLE_TOO_LONG:
          case State::ASSEMBLE_SUCCESS:
            return forward_read;
          case State::ASSEMBLE_TOO_SHORT:
            if (forward_read.size() > longest_read.size()) {
              longest_read = std::move(forward_read);
            }
            break;
          case State::ASSEMBLE_FAILED:
            break;
        }

        auto reverse_read = graph.find_read<true>(sink, source);
        switch (get_read_state(reverse_read)) {
          case State::ASSEMBLE_TOO_LONG:
          case State::ASSEMBLE_SUCCESS:
            return reverse_read;
          case State::ASSEMBLE_TOO_SHORT:
            if (reverse_read.size() > longest_read.size()) {
              longest_read = std::move(reverse_read);
            }
            break;
          case State::ASSEMBLE_FAILED:
            break;
        }
      }
    }
    return longest_read;
  }

 public:
  ReadAssembler(const std::size_t read_len) { this->read_len = read_len; }

  // delete copy constructor and copy assignment
  ReadAssembler(const ReadAssembler&) = delete;

  ReadAssembler& operator=(const ReadAssembler&) = delete;

  /* move constructor */
  ReadAssembler(ReadAssembler&& rhs) {
    read_len = rhs.read_len;
    sequences = std::move(rhs.sequences);
  }

  /* move assignment */
  ReadAssembler& operator=(ReadAssembler&& rhs) {
    read_len = rhs.read_len;
    sequences = std::move(rhs.sequences);
    return *this;
  }

  /**
   * @brief Adds a sequence to the list of sequences.
   *
   * This function adds a sequence to the list of sequences. If the type of the
   * sequence is bio::istring, it is directly added to the list. Otherwise, a
   * new Sequence<bio::istring> object is created with the provided sequence
   * information and added to the list.
   *
   * @tparam T The type of the sequence.
   * @param seq The sequence to be added.
   */
  template <class T>
  auto add_seq(const Sequence<T>& seq) {
    if constexpr (std::is_same_v<T, bio::istring>) {
      sequences.push_back(seq);
    } else {
      sequences.emplace_back(Sequence<bio::istring>{
          seq.read_id, seq.left_bound, seq.right_bound,
          bio::Codec::to_istring(seq.seq), seq.qual, seq.forward_strain});
    }
  }

  auto assemble() {
    auto kmer_size = std::ceil(std::log2(read_len) * 2.5);

    /* For each iteration, increase kmer_size */
    for (std::size_t i = 0; i < param.MAX_ITERATIONS; ++i) {
      // spdlog::debug("Assemble: kmer_size = {}", kmer_size);
      auto failed = false;
      auto corrected_read = assemble_read(kmer_size);

      switch (get_read_state(corrected_read)) {
        case State::ASSEMBLE_SUCCESS:
          return corrected_read;
        case State::ASSEMBLE_TOO_LONG:
        case State::ASSEMBLE_TOO_SHORT:
          // spdlog::debug("Assemble read length({}) is not ideal({} - {}),
          // increase "
          //               "kmer size({})",
          //               corrected_read.size(),
          //               std::size_t(read_len *
          //               param.ASSEMBLE_TOO_SHORT_RATIO), std::size_t(read_len
          //               * param.ASSEMBLE_MAX_LEN_RATIO), kmer_size +
          //               param.INCREASE_KMER_SIZE);
          kmer_size += param.INCREASE_KMER_SIZE;
          break;
        case State::ASSEMBLE_FAILED:
          failed = true;
          // spdlog::debug("Assemble failed");
          break;
      }
      if (failed) {
        break;
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

  auto clear() { sequences.clear(); }

 private:
  std::size_t read_len;
  std::vector<Sequence<bio::istring>> sequences;
};