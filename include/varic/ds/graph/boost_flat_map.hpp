#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <numeric>
#include <ranges>
#include <set>
#include <span>
#include <stack>
#include <string>
#include <string_view>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <biovoltron/algo/assemble/graph/graph_wrapper.hpp>

namespace bio = biovoltron;
namespace fs = std::filesystem;

using ByteUnit = std::uint64_t;

constexpr static auto BASES_PER_BYTE = 4u;
constexpr static auto BASES_PER_BYTE_UNIT = BASES_PER_BYTE * sizeof(ByteUnit);
constexpr static auto MAX_SPAN_LENGTH = 2u;
constexpr static auto MAX_KMER_SIZE = BASES_PER_BYTE_UNIT * MAX_SPAN_LENGTH;

/* three-way operator for std::span<> */
// template <typename T>
// auto operator<=>(const std::span<T>& a, const std::span<T>& b) {
//   return a.data() <=> b.data();
// }

// TODO: reallocate cause the memory address broken, need to fix it
struct TestGraph {
 public:
  struct Param {
    /* For global assmebling */
    /* The valid weight ratio that should consider as a valid branch */
    const double VALID_BRANCH_RATIO = 0.2;

    /* For local assembling */
    /* The weight ratio of a branch point to the other */
    const double LOCAL_ASSEMBLE_WEIGHT_RATIO = 0.8;
    /* The length threshold for preventing the cycle case */
    const double LOCAL_ASSEMBLE_LENGTH_RATIO = 1.2;
  } param;

  // use std::array or std::bitset?
  using Kmer = std::array<ByteUnit, MAX_SPAN_LENGTH>;
  using KmerView = std::span<ByteUnit>;

  struct KmerViewHash {
    std::size_t operator()(const KmerView& kmer) const {
      return std::accumulate(kmer.begin(), kmer.end(), ByteUnit(0),
                            [](const auto& a, const auto& b) {
                              return a ^ std::hash<ByteUnit>{}(b);
                            });
    }
  };

  struct KmerViewEqual {
    bool operator()(const KmerView& a, const KmerView& b) const {
      return std::ranges::equal(a, b);
    }
  };

  struct KmerViewLess {
    bool operator()(const KmerView& a, const KmerView& b) const {
      return std::ranges::lexicographical_compare(a, b);
    }
  };

  struct KmerLess {
    bool operator()(const Kmer& a, const Kmer& b) const {
      return std::ranges::lexicographical_compare(a, b);
    }
  };

  // if you use your own flat_map, you need to maintain the KmerView 
  template<class T>
  // using KmerMap = std::unordered_map<KmerView, T, KmerViewHash, KmerViewEqual>;
  using KmerMap = boost::container::flat_map<Kmer, T, KmerLess>;


  struct VertexProperty {
    // index at kmers
    std::size_t index;
    std::set<std::size_t> appearances;
  };

  struct EdgeProperty {
    std::size_t count = 0;
  };

  using Graph = typename bio::GraphWrapper<VertexProperty, EdgeProperty>;
  using Vertex = typename Graph::Vertex;
  using Edge = typename Graph::Edge;
  using Path = std::vector<std::pair<Vertex, std::size_t>>;

  auto create_vertex(const Kmer& kmer) {

    const auto v = g.create_vertex();
    kmer_to_vertex[kmer] = v;

    // debug
    // vertex_to_kmer[v] = kmer;
    return v;
  }

  auto create_edge(const Vertex& u, const Vertex& v) {
    const auto e = g.create_edge(u, v);
    // debug
    edge_to_weight[e] = g[e].count;
    return e;
  }

  auto get_vertex(const Kmer& kmer) {
    if (kmer_to_vertex.contains(kmer)) {
      return kmer_to_vertex[kmer];
    }
    return create_vertex(kmer);
  }

  auto get_edge(const Vertex& u, const Vertex& v) {
    for (const auto e : g.out_edges(u, false)) {
      if (g.target(e) == v) {
        return e;
      }
    }
    return create_edge(u, v);
  }

  auto increase_edge_weight(const Vertex& u, const Vertex& v,
                            std::size_t weight = 1) {
    const auto e = get_edge(u, v);
    g[e].count += weight;
  }

  auto out_weight_sum(const Vertex& v) {
    auto out_edges = g.out_edges(v, false);
    return std::accumulate(
        out_edges.begin(), out_edges.end(), 0.0,
        [this](const auto& a, const auto& b) { return a + g[b].count; });
  }

  auto in_weight_sum(const Vertex& v) {
    auto in_edges = g.in_edges(v, false);
    return std::accumulate(
        in_edges.begin(), in_edges.end(), 0.0,
        [this](const auto& a, const auto& b) { return a + g[b].count; });
  }

  auto get_vertex_inside_range(const std::size_t start, const std::size_t end) {
    assert(start < end);
    auto vertices = std::vector<Vertex>{};
    for (const auto& [kmer, v] : kmer_to_vertex) {
      const auto& appearances = g[v].appearances;
      auto it = appearances.lower_bound(start);
      if (it != appearances.end() && *it <= end) {
        vertices.push_back(v);
      }
    }
    // for (const auto& [kmer, appearance] : kmer_appearances) {
    //   auto it = appearance.lower_bound(start);
    //   if (it != appearance.end() && *it <= end) {
    //     vertices.push_back(kmer_to_vertex[kmer]);
    //   }
    // }
    return vertices;
  }

  auto is_source(const Vertex& v, const std::size_t weight_threshold) {
    auto out_degree = std::ranges::count_if(
        g.out_edges(v, false),
        [&, this](const auto& e) { return g[e].count >= weight_threshold; });
    auto in_degree = std::ranges::count_if(
        g.in_edges(v, false),
        [&, this](const auto& e) { return g[e].count >= weight_threshold; });
    return in_degree == 0 && out_degree > 0;
  }

  auto is_source(const Kmer& kmer, const std::size_t weight_threshold) {
    auto v = get_vertex(kmer);
    return is_source(v, weight_threshold);
  }

  auto is_sink(const Vertex& v, const std::size_t weight_threshold) {
    auto out_degree = std::ranges::count_if(
        g.out_edges(v, false),
        [&, this](const auto& e) { return g[e].count >= weight_threshold; });
    auto in_degree = std::ranges::count_if(
        g.in_edges(v, false),
        [&, this](const auto& e) { return g[e].count >= weight_threshold; });
    return in_degree > 0 && out_degree == 0;
  }

  auto is_sink(const Kmer kmer, const std::size_t weight_threshold) {
    auto v = get_vertex(kmer);
    return is_sink(v, weight_threshold);
  }

  /**
   * @brief Call this when the graph is ready to assemble something
   * 
   * @return auto 
   */
  auto set_vertex_index() {
    for (auto idx = 0; auto [kmer, v] : kmer_to_vertex) {
      assert(idx == kmer_to_vertex.index_of(kmer_to_vertex.find(kmer)));
      g[v].index = idx;
      idx += 1;
    }
  }

  auto get_kmer(const Vertex& v) -> KmerView {
    return kmer_to_vertex.nth(g[v].index)->first;
    // return kmers[g[v].index];
  }

  auto add_pos_to_vertex(const Vertex& v, const std::size_t pos) {
    g[v].appearances.insert(pos);
  }

  template <bool Reverse = false>
  auto concat_vertices(const Vertex& source, Path& path) {
    assert(valid_path<Reverse>(path));
    auto read = decode_kmer_to_seq(get_kmer(source));
    for (const auto& [v, weight] : path) {
      auto now_kmer = get_kmer(v);
      if constexpr (Reverse == false) {
        read += get_last_base(now_kmer);
      } else {
        read += get_first_base(now_kmer);
      }
    }
    if constexpr (Reverse) {
      std::ranges::reverse(read);
    }
    return bio::Codec::to_string(read);
  }

  // TODO: extend source and sink
  template <bool Reverse = false>
  auto extend_path_at_end(Path& path) {}

  /* do local assemble between source and sink */
  template <bool Reverse = false>
  auto local_assemble(const Vertex& source, const Vertex& sink,
                      const std::size_t length_limit) {
    // spdlog::debug("Local assemble between {} - {}", g[source].kmer,
    //               g[sink].kmer);
    // TODO: considering following stuffs
    // 1. cycle
    // 2. the possibility of jumping to vertexes behind the sink
    //   - set a threshold of the path length
    auto paths = std::vector<Path>{};
    auto dfs = [&](auto&& self, const Vertex& now, Path& path,
                   std::set<Vertex>& vis) {
      if (now == sink) {
        paths.push_back(path);
        return;
      }
      if (path.size() >= length_limit) {
        return;
      }
      vis.insert(now);
      auto edges = std::vector<Edge>{};
      if constexpr (Reverse == false) {
        edges = g.out_edges(now, false);
      } else {
        edges = g.in_edges(now, false);
      }
      for (const auto& e : edges) {
        auto u = Vertex{};
        if constexpr (Reverse == false) {
          u = g.target(e);
        } else {
          u = g.source(e);
        }
        if (vis.contains(u)) {
          continue;
        }
        path.emplace_back(u, g[e].count);
        self(self, u, path, vis);
        path.pop_back();
      }
      vis.erase(now);
    };
    /* the path doesn't contain source for preserving weight between source and
     * its previous vertex */
    auto path = Path{};
    auto vis = std::set<Vertex>{};
    dfs(dfs, source, path, vis);
    // sort the paths by average weight
    std::ranges::sort(paths, [](const auto& a, const auto& b) {
      auto a_sum = std::accumulate(
          a.begin(), a.end(), 0.0,
          [](const auto& a, const auto& b) { return a + b.second; });
      auto b_sum = std::accumulate(
          b.begin(), b.end(), 0.0,
          [](const auto& a, const auto& b) { return a + b.second; });
      return a_sum / a.size() > b_sum / b.size();
    });
    // spdlog::debug("There are {} paths.", paths.size());
    // for (auto& path : paths) {
      // auto weight_sum = std::accumulate(
      //     path.begin(), path.end(), 0.0,
      //     [](const auto &a, const auto &b) { return a + b.second; });
      // auto [mn, mx] = std::ranges::minmax_element(
      //     path.begin(), path.end(),
      //     [](const auto &a, const auto &b) { return a.second < b.second; });
      // auto s = std::string{};
      // for (auto [v, w] : path) {
      //   s.push_back(g[v].kmer.back());
      // }
      // spdlog::debug("candidate({}, {}, {:.3f}): {}, ", mn->second,
      // mx->second,
      //               weight_sum / path.size(), s);
    // }
    assert(paths.size() > 0);
    return paths[0];
  }

  template <bool Reverse = false>
  auto valid_path(const Path& path) {
    auto sz = path.size();

    auto print_path = [&]() {
      for (int i = 0; i < path.size(); i++) {
        auto [v, weight] = path[i];
        auto kmer = get_kmer(v);
        spdlog::debug("i = {}, v = {}, data = {}, weight = {}", i,
                      bio::Codec::to_string(decode_kmer_to_seq(kmer)), fmt::ptr(kmer.data()), weight);
      }
    };

    for (auto i = 1u; i < sz; i++) {
      auto [prev, prev_weight] = path[i - 1];
      auto [now, now_weight] = path[i];
      auto prev_kmer = get_kmer(prev);
      auto now_kmer = get_kmer(now);
      if constexpr (Reverse == false) {
        auto a = decode_kmer_to_seq(prev_kmer).substr(1);
        auto b = decode_kmer_to_seq(now_kmer).substr(0, kmer_size - 1);
        if (a != b) {
          spdlog::debug("Rev = {}, i = {}, prev = {}, now = {}", Reverse, i,
                        bio::Codec::to_string(a), 
                        bio::Codec::to_string(b));
          print_path();
          return false;
        }
      } else {
        auto a = decode_kmer_to_seq(prev_kmer).substr(0, kmer_size - 1);
        auto b = decode_kmer_to_seq(now_kmer).substr(1);
        if (a != b) {
          spdlog::debug("Rev = {}, i = {}, prev = {}, now = {}", Reverse, i,
                        bio::Codec::to_string(a), bio::Codec::to_string(b));
          print_path();
          return false;
        }
      }
    }
    return true;
  }

  // check the path for the need of local reassemble
  template <bool Reverse = false>
  auto check_and_reassemble_path(const Path& path) {
    // 1. identity the source and sink that need local assembly
    //   - if the adjancent vertexes in the path are not supported by each
    //   other,
    //   - then it could need local reassmebly

    // for (auto [v, w] : path) {
    //   spdlog::debug("v = {}, weight = {}", g[v].kmer, w);
    // }

    auto new_path = Path{};
    new_path.reserve(path.size());
    auto find_previous_breakpoint = [&](const std::size_t weight) {
      for (auto i = (int) new_path.size() - 1; i >= 0; i--) {
        if (new_path[i].second > weight) {
          return i;
        }
      }
      return -1;
    };
    assert(valid_path<Reverse>(path));

    for (auto i = 0u; i < path.size(); i++) {
      if (i && path[i].second > 2 * path[i - 1].second) {
        // spdlog::debug("Found branch point at {}, i = {}",
        // g[path[i].first].kmer, i);
        auto prev_branch_pos = find_previous_breakpoint(
            path[i].second * param.LOCAL_ASSEMBLE_WEIGHT_RATIO);
        if (prev_branch_pos != -1) {
          auto prev_branch_v = new_path[prev_branch_pos].first;
          auto old_path_len = new_path.size() - prev_branch_pos;
          if (old_path_len >= 2 * kmer_size) {
            // the path is too long, so the reassembled process may take a lot
            // time, skip it
            new_path.emplace_back(path[i]);
            continue;
          }
          auto now_branch_v = path[i].first;
          // spdlog::debug("Found previous branch point at {}, len = {}",
          //               g[prev_branch_v].kmer, old_path_len);
          assert(valid_path<Reverse>(new_path));
          // spdlog::debug("Found source point at {}, len = {}",
          //               g[prev_branch_v].kmer, old_path_len);
          auto reassembled_path = local_assemble<Reverse>(
              prev_branch_v, now_branch_v,
              old_path_len * param.LOCAL_ASSEMBLE_LENGTH_RATIO);
          // spdlog::debug("Found sink point at {}, len = {}",
          // g[now_branch_v].kmer,
          //               reassembled_path.size());
          // spdlog::debug("=========================\n\n");

          /* Remove old path, add reassembled path */
          while (new_path.size() != prev_branch_pos + 1) {
            // spdlog::debug("pop {}", g[new_path.back().first].kmer);
            new_path.pop_back();
          }
          std::ranges::copy(reassembled_path, std::back_inserter(new_path));
          // for (auto [v, w] : reassembled_path) {
          //   spdlog::debug("Re: v = {}, weight = {}", g[v].kmer, w);
          // }
          assert(valid_path<Reverse>(new_path));
        } else {
          // spdlog::debug("New source found, need extend it!");
          new_path.clear();
          new_path.emplace_back(path[i]);
          // extend from source
        }
      } else {
        new_path.emplace_back(path[i]);
        // spdlog::debug("push {}", g[path[i].first].kmer);
        // spdlog::debug("new_path:");
        // for (auto [v, w] : new_path | std::views::reverse |
        // std::views::take(10) | std::views::reverse) {
        //   spdlog::debug("v = {}, weight = {}", g[v].kmer, w);
        // }
        // assert(valid_path<Reverse>(new_path));
      }
    }
    assert(valid_path<Reverse>(new_path));
    return new_path;

    // TODO: remove low weight path at the end
  }

  template <bool Reverse = false>
  auto path_finder(const Vertex& source, const Vertex& sink, Path& path) {
    auto vis = std::map<Vertex, std::size_t>{};
    auto fail_vertexes = std::set<Vertex>{};

    /* record (Vertex, edge_weight, index of vertex in path) */
    auto stk = std::stack<std::tuple<Vertex, std::size_t, std::size_t>>{};

    path.reserve(max_assemble_len);
    stk.emplace(source, 0, 0);

    int branch_cnt = 0;
    while (!stk.empty()) {
      const auto [v, w, pos] = stk.top();
      stk.pop();
      if (v == sink) {
        // spdlog::debug("Reach sink {}, branch_cnt = {}", g[v].kmer,
        // branch_cnt);
        return true;
      }
      if (pos >= max_assemble_len) {
        // spdlog::debug("Reach max_assemble_len {}", max_assemble_len);
        continue;
      }
      if (fail_vertexes.contains(v)) {
        continue;
      }
      bool need_pop = false;
      while (path.size() > pos) {
        if (path.size() > max_path_size) {
          max_path_size = path.size();
          /* record the path */
          max_path = path;
        }
        auto [u, w] = path.back();
        fail_vertexes.insert(u);
        path.pop_back();
        need_pop = true;
      }

      path.emplace_back(v, w);
      if (need_pop) {
        assert(valid_path<Reverse>(path));
      }

      vis[v] += 1;
      auto edges = std::vector<Edge>{};
      if constexpr (Reverse == false) {
        edges = g.out_edges(v, false);
      } else {
        edges = g.in_edges(v, false);
      }
      auto edge_weights = std::vector<std::size_t>(edges.size());
      auto edge_weight_sum = 0.0;
      for (auto i = 0u; const auto& e : edges) {
        auto u = Vertex{};
        if constexpr (Reverse == false) {
          u = g.target(e);
        } else {
          u = g.source(e);
        }
        // edge_weights[i] = g[e].count / std::max(1ul, vis[u]);
        edge_weights[i] = g[e].count / (vis[u] + 1);
        edge_weight_sum += edge_weights[i];
        i += 1;
      }

      auto order = std::vector<std::size_t>(edges.size());
      std::iota(order.begin(), order.end(), 0);
      std::ranges::sort(order, [&](const auto& a, const auto& b) {
        return edge_weights[a] > edge_weights[b];
      });

      auto cnt = 0;
      for (const auto x : order) {
        auto w = edge_weights[x];
        if (w == 0) {
          continue;
        }
        if (w >= param.VALID_BRANCH_RATIO * edge_weight_sum) {
          cnt += 1;
        }
      }
      if (cnt != 1) {
        branch_cnt += 1;
      }

      /* we wanna let the heaviest edge at the top of stack, so reverse order*/
      for (const auto x : order | std::views::reverse) {
        auto [e, w] = std::tie(edges[x], edge_weights[x]);
        auto u = Vertex{};
        if constexpr (Reverse == false) {
          u = g.target(e);
        } else {
          u = g.source(e);
        }
        // w == 0 -> vis[u] > weight -> impossible
        if (w == 0 || w < param.VALID_BRANCH_RATIO * edge_weight_sum) {
          continue;
        }
        stk.emplace(u, w, pos + 1);
      }
    };
    return false;
  };

  auto reset_vertex_weight_on_path(Path& path) {
    // set the duplicated kmer weight to the minimum weight it appear
    auto min_weight = std::map<Vertex, std::size_t>{};
    for (const auto& [vertex, weight] : path) {
      if (min_weight.contains(vertex)) {
        min_weight[vertex] = std::min(min_weight[vertex], weight);
      } else {
        min_weight[vertex] = weight;
      }
    }
    for (auto& [vertex, weight] : path) {
      weight = min_weight[vertex];
    }
  }

 public:
  /**
   * @brief Constructs a TestGraph object with the given parameters.
   *
   * @param kmer_size The size of the k-mer.
   * @param max_assemble_len The maximum length for assembly.
   */
  TestGraph(const std::size_t kmer_size, const std::size_t max_assemble_len) {
    if (kmer_size > MAX_KMER_SIZE) {
      spdlog::warn("kmer_size should be less than {}", MAX_KMER_SIZE);
      this->kmer_size = MAX_KMER_SIZE;
    } else {
      this->kmer_size = kmer_size;
    }
    this->used_bytes = (kmer_size + BASES_PER_BYTE_UNIT - 1) / BASES_PER_BYTE_UNIT;
    this->shift_at_last_unit = (kmer_size - 1) % BASES_PER_BYTE_UNIT * 2;
    this->max_assemble_len = max_assemble_len;
    this->kmer_to_vertex.reserve(this->max_assemble_len * 2);
  }

  auto get_kmer_size() const { return kmer_size; }

  /* kmer related function */

  auto pretty_print(const KmerView kmer) {
    auto seq = decode_kmer_to_seq(kmer);
    return bio::Codec::to_string(seq);
  }

  auto pretty_print(const ByteUnit unit) {
    std::stringstream ss;
    ss << std::hex << unit;
    return ss.str();
  }

  auto decode_unit_to_seq(const ByteUnit unit) -> bio::istring {
    auto seq = bio::istring{};
    for (auto i = 0u; i < BASES_PER_BYTE_UNIT; i++) {
      seq.push_back((unit >> (i * 2)) & 3);
    }
    return seq;
  }

  auto decode_kmer_to_seq(const KmerView kmer) -> bio::istring {
    auto seq = bio::istring{};
    for (auto i = 0u; i < used_bytes; i++) {
      seq += decode_unit_to_seq(kmer[i]);
    }
    seq = seq.substr(0, kmer_size);
    return seq;
  }

  auto encode_seq_to_unit(const bio::istring_view seq) {
    assert(seq.size() <= BASES_PER_BYTE_UNIT);
    auto unit = ByteUnit{0};
    for (auto i = 0u; i < std::min(seq.size(), BASES_PER_BYTE_UNIT); i++) {
      unit |= ByteUnit(seq[i]) << (i * 2);
    }
    return unit;
  }

  auto encode_seq_to_kmer(const bio::istring_view seq) {
    assert(seq.size() == kmer_size);
    auto kmer = Kmer{};
    for (auto i = 0u; i < kmer_size; i += BASES_PER_BYTE_UNIT) {
      auto subseq = seq.substr(i, BASES_PER_BYTE_UNIT);
      auto data = encode_seq_to_unit(seq.substr(i, BASES_PER_BYTE_UNIT));
      // spdlog::debug("encode {} to {}", bio::Codec::to_string(subseq), pretty_print(data));
      kmer[i / BASES_PER_BYTE_UNIT] = data;
    }
    // spdlog::debug("kmer = {}", pretty_print(kmer));
    return kmer;
  }

  auto get_base_at_unit(const ByteUnit unit, const std::size_t pos) {
    assert(pos < BASES_PER_BYTE_UNIT);
    return (unit >> (pos * 2)) & 3;
  }

  auto get_first_base(const KmerView kmer) -> bio::ichar {
    return kmer[0] & 0b11;
  }

  auto get_last_base(const KmerView kmer) -> bio::ichar {
    return kmer[used_bytes - 1] >> shift_at_last_unit;
  }

  /**
   * @brief Shifts the given Kmer and appends the specified base.
   * 
   * @param kmer The Kmer to be shifted and appended.
   * @param base The base to be appended to the Kmer.
   * @return The updated Kmer after shifting and appending the base.
   */
  auto shift_and_append(Kmer kmer, const ByteUnit base) {
    for (auto i = 0u; i < used_bytes; i++) {
      kmer[i] = kmer[i] >> 2;
      if (i + 1 < used_bytes) {
        kmer[i] |= (kmer[i + 1] & 0b11) << ((BASES_PER_BYTE_UNIT - 1) * 2);
      } else {
        kmer[i] |= base << shift_at_last_unit;
      } 
    }
    return kmer;
  }

  auto add_seq(bio::istring_view sequence, const std::size_t left_bound) {
    if (sequence.size() < kmer_size) {
      return;
    }

    auto prev_kmer = encode_seq_to_kmer(sequence.substr(0, kmer_size));
    auto prev_vertex = get_vertex(prev_kmer);


    // don't use `prev_kmer`! we use Kmer.data() for comparing
    // if we use `prev_kmer` directly, then we need to change compare criteria
    // set_kmer_pos(g[prev_kmer].kmer, 0);
    add_pos_to_vertex(prev_vertex, 0 + left_bound);
    
    for (auto i = 1u; i + kmer_size <= sequence.size(); i++) {
      auto now_kmer = shift_and_append(prev_kmer, sequence[i + kmer_size - 1]);
      assert(decode_kmer_to_seq(now_kmer) == sequence.substr(i, kmer_size));
      auto now_vertex = get_vertex(now_kmer);
      add_pos_to_vertex(now_vertex, i + left_bound);
      increase_edge_weight(prev_vertex, now_vertex);
      prev_vertex = now_vertex;
      prev_kmer = now_kmer;
    }
  }

  
  auto is_source(bio::istring_view seq, const std::size_t weight_threshold) {
    if (seq.size() != kmer_size) {
      spdlog::warn("The size of the sequence should be {}, not {}", kmer_size, seq.size());
      return false;
    }
    auto kmer = encode_seq_to_kmer(seq);
    return is_source(kmer, weight_threshold);
  }

  auto is_sink(bio::istring_view seq, const std::size_t weight_threshold) {
    if (seq.size() != kmer_size) {
      spdlog::warn("The size of the sequence should be {}, not {}", kmer_size, seq.size());
      return false;
    }
    auto kmer = encode_seq_to_kmer(seq);
    return is_sink(kmer, weight_threshold);
  }

  /**
   * @brief add one sequence to the graph
   * @param sequence the sequence to be added
   * @param left_bound the left boundary of the sequence in the read
   */
  // auto add_seq(const std::string& sequence, const std::size_t left_bound) {
  //   if (sequence.size() < kmer_size) {
  //     return;
  //   }

  //   auto update_appearance = [&](KmerView kmer, const std::size_t pos) {
  //     kmer_appearances[kmer].insert(pos + left_bound);
  //   };

  //   auto prev_vertex = get_vertex(sequence.substr(0, kmer_size));
  //   update_appearance(prev_vertex, 0);
  //   for (auto i = 1u; i + kmer_size <= sequence.size(); i++) {
  //     bool found = false;
  //     for (const auto& edge : g.out_edges(prev_vertex, false)) {
  //       const auto u = g.target(edge);
  //       if (g[u].kmer.back() == sequence[i + kmer_size - 1]) {
  //         g[edge].count += 1;
  //         update_appearance(u, i);
  //         prev_vertex = u;
  //         found = true;
  //         break;
  //       }
  //     }
  //     if (!found) {
  //       const auto v = get_vertex(sequence.substr(i, kmer_size));
  //       update_appearance(v, i);
  //       create_edge(prev_vertex, v);
  //       prev_vertex = v;
  //     }
  //   }
  // }

  auto get_sources(const std::size_t start, const std::size_t end,
                   const std::size_t weight_threshold) {
    auto candidate_sources = get_vertex_inside_range(start, end);
    // spdlog::debug("Find Sources({} - {}) using weight = {}: candidate.size() = {}", start, end, weight_threshold,
    //               candidate_sources.size());
    auto sources = std::vector<Vertex>{};
    std::ranges::copy_if(candidate_sources, std::back_inserter(sources),
                         [this, weight_threshold](const auto& v) {
      return is_source(v, weight_threshold);
    });
    std::ranges::sort(sources, [this](const auto& a, const auto& b) {
      const auto& a_appearances = g[a].appearances;
      const auto& b_appearances = g[b].appearances;
      assert(a_appearances.size() > 0 && b_appearances.size() > 0);
      
      const auto a_first_pos = *a_appearances.begin();
      const auto a_last_pos = *a_appearances.rbegin();
      const auto b_first_pos = *b_appearances.begin();
      const auto b_last_pos = *b_appearances.rbegin();
      
      const auto a_diff = a_last_pos - a_first_pos;
      const auto b_diff = b_last_pos - b_first_pos;
      if (a_diff != b_diff) {
        return a_diff < b_diff;
      }
      if (a_last_pos != b_last_pos) {
        return a_last_pos < b_last_pos;
      }
      return b_first_pos > a_first_pos;
    });
    return sources;
  }

  auto get_sinks(const std::size_t start, const std::size_t end,
                 const std::size_t weight_threshold) {
    auto candidate_sinks = get_vertex_inside_range(start, end);
    // spdlog::debug("Find Sinks({} - {}) using weight = {}: candidate.size() = {}", start, end, weight_threshold,
    //               candidate_sinks.size());
    auto sinks = std::vector<Vertex>{};
    std::ranges::copy_if(candidate_sinks, std::back_inserter(sinks),
                         [this, weight_threshold](const auto& v) {
      return is_sink(v, weight_threshold);
    });
    std::ranges::sort(sinks, [this](const auto& a, const auto& b) {
      const auto& a_appearances = g[a].appearances;
      const auto& b_appearances = g[b].appearances;
      assert(a_appearances.size() > 0 && b_appearances.size() > 0);

      const auto a_first_pos = *a_appearances.begin();
      const auto a_last_pos = *a_appearances.rbegin();
      const auto b_first_pos = *b_appearances.begin();
      const auto b_last_pos = *b_appearances.rbegin();

      const auto a_diff = a_last_pos - a_first_pos;
      const auto b_diff = b_last_pos - b_first_pos;
      if (a_diff != b_diff) {
        return a_diff < b_diff;
      }
      if (a_first_pos != b_first_pos) {
        return a_first_pos > b_first_pos;
      }
      return a_last_pos < b_last_pos;
    });
    return sinks;
  }

  /**
   * @brief Get the read object
   *
   * @param source
   * @param sink ideal sink
   * @return auto
   */
  template <bool Reverse = false>
  auto find_read(const Vertex& source, const Vertex& sink) {
    assert(source != sink);
    // spdlog::debug("Find the path in {} direction", Reverse ? "reverse" :
    // "forward"); spdlog::debug("Try source = {}({} - {}) sink = {}({} - {})",
    // g[source].kmer,
    //               *g[source].appearances.begin(),
    //               *g[source].appearances.rbegin(), g[sink].kmer,
    //               *g[sink].appearances.begin(),
    //               *g[sink].appearances.rbegin());

    set_vertex_index();

    /* reset max_path */
    max_path_size = 0;
    max_path.clear();

    /* store (Vertex, edge_weight) for further locally reassemble */
    auto path = Path{};

    // auto vis = std::map<Vertex, std::size_t>{};
    // auto fail_vertexes = std::set<Vertex>{};
    // auto result = path_finder<Reverse>(source, sink, path, vis,
    // fail_vertexes); if (result == false) {
    //   if (max_path_size == max_assemble_len) {
    //     // spdlog::debug("max_path length reach threshold {}",
    //     max_assemble_len); return concat_vertices<Reverse>(source, max_path);
    //   }
    //   path = max_path;
    // }

    auto result = path_finder<Reverse>(source, sink, path);
    if (result == false) {
      // spdlog::debug("Cannot find a proper path for given source and sink");
      // spdlog::debug("Using path with maximum length({}) as result",
      //               max_path_size);
      /* if stuck in cycle, need increase kmer_size, return directly */
      if (max_path_size == max_assemble_len) {
        // spdlog::debug("max_path length reach threshold {}",
        // max_assemble_len);
        return concat_vertices<Reverse>(source, max_path);
      }

      // TODO: trim the path
      // Since we didn't reach the sink, the end of path should have less
      // confidence, need to do something like trimmng
      // 1. find the last vertex that has weight > 1
      //   - may reassemble from last vertex
      //

      // for (auto [v, w] : max_path | std::views::reverse |
      // std::views::take(100) | std::views::reverse) {
      //   spdlog::debug("max: v = {}, w = {}", g[v].kmer, w);
      // }
      path = max_path;
    }

    // spdlog::debug("before path.size() = {}", path.size());
    // assert(valid_path<Reverse>(path));
    reset_vertex_weight_on_path(path);
    // assert(valid_path<Reverse>(path));

    path = check_and_reassemble_path<Reverse>(path);
    // print_graphviz(path);
    // spdlog::debug("after path.size() = {}", path.size());
    return concat_vertices<Reverse>(source, path);

    // auto fail_vertexes = std::set<Vertex>{};
    // auto stk = std::stack<Vertex>{};

    // stk.emplace(source);

    // while (1) {
    //   if (stk.empty()) {
    //     break;
    //   }
    //   auto v = stk.top();
    //   stk.pop();
    //   // {
    //   //   spdlog::debug("v = {}, ({} - {}), dep = {}", g[v].kmer,
    //   *g[v].appearances.begin(),
    //   //                 *g[v].appearances.rbegin(), path.size());
    //   //   // spdlog::debug("cnt = {}", cnt);
    //   // }
    //   path.emplace_back(v);
    //   if (v == sink) {
    //     spdlog::debug("Reach sink {}", g[v].kmer);
    //     break;
    //   }
    //   /* if path.size() >= max_len -> stuck in cycle */
    //   if (path.size() >= max_len) {
    //     spdlog::debug("Reach max_len {}", max_len);
    //     break;
    //   }
    //   vis[v] += 1;
    //   auto out_edges = g.out_edges(v, false);
    //   std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b)
    //   {
    //     const auto a_vis_times = std::max(1, vis[g.target(a)]);
    //     const auto b_vis_times = std::max(1, vis[g.target(b)]);
    //     return static_cast<double>(g[a].count) / a_vis_times >
    //            static_cast<double>(g[b].count) / b_vis_times;
    //   });
    //   for (const auto& e : out_edges) {
    //     auto u = g.target(e);
    //     if (vis[u] >= g[e].count) {
    //       continue;
    //     }
    //     // Don't use this, will skip some vertices when the region is complex
    //     // if (g[e].count == 1) {
    //     //   continue;
    //     // }
    //     spdlog::debug("now = {}, g[e].kmer count = {}", g[v].kmer,
    //     g[e].count); stk.push(u); break;
    //   }
    // }

    // TODO: we may trim the path by the weight of edge
    // 1. modify the path structure for storing the edge
    // 2. For local reassemble, You need to do topologically sort first, get
    //    the order, then parse the path between two nodes(source and sink)
    //    locally.
    // for (auto last = path[0]; const auto& v : path | std::views::drop(1)) {
    //   if (v == last) {
    //     continue;
    //   }
    //   auto out_edges = g.out_edges(last, false);
    //   auto weight_sum = std::accumulate(out_edges.begin(), out_edges.end(),
    //   0.0,
    //                                     [this](const auto &a, const auto &b)
    //                                     {
    //                                       return a + g[b].count;
    //                                     });
    //   auto weight_ratio = std::vector<double>(out_edges.size());
    //   auto branch_cnt = 0;
    //   for (auto i = 0u; const auto& e : g.out_edges(last, false)) {
    //     weight_ratio[i++] = g[e].count / weight_sum;
    //     if (weight_ratio >= 0.5) {

    //     }
    //   }
    //   last = v;
    // }

    // return concat_vertices(path);

    // auto stk = std::stack<std::pair<Vertex, std::size_t>>{};
    // stk.emplace(source, *g[source].appearances.begin());
    // while (!stk.empty()) {
    //   const auto [v, pos] = stk.top();
    //   stk.pop();

    //   path.push_back(v);
    //   if (v == sink) {
    //     spdlog::debug("Reach sink {}", g[v].kmer);
    //     break;
    //   }
    //   vis[v] += 1;
    //   auto out_edges = g.out_edges(v, false);
    //   // TODO: can be imporved by use vector of edges -> query vis only once
    //   std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b)
    //   {
    //     const auto a_vis_times = std::max(1, vis[g.target(a)]);
    //     const auto b_vis_times = std::max(1, vis[g.target(b)]);
    //     return g[a].count / a_vis_times > g[b].count / b_vis_times;
    //   });

    //   // remember the visis times and divide by visit time?
    //   // take the
    //   auto is_nearby = [&](const std::set<std::size_t> &appearances) {
    //     auto supposed_pos = pos + 1;
    //     auto it = appearances.lower_bound(supposed_pos);
    //     if (it == appearances.begin()) {
    //       if (*it - supposed_pos > param.NEARBY_THRESHOLD) {
    //         return static_cast<std::size_t>(-1);
    //       }
    //       // return *it;
    //       return (*it + supposed_pos) / 2;
    //     }

    //     // auto prev_it = std::prev(it);
    //     // if (it == appearances.end()) {
    //     //   if (supposed_pos - *prev_it > kmer_size) {
    //     //     return static_cast<std::size_t>(-1);
    //     //   }
    //     //   return *prev_it;
    //     // }

    //     // if (supposed_pos - *prev_it > kmer_size &&
    //     //     *it - supposed_pos > kmer_size) {
    //     //   return static_cast<std::size_t>(-1);
    //     // }

    //     // if (supposed_pos - *prev_it < *it - supposed_pos) {
    //     //   return (*prev_it + supposed_pos) / 2;
    //     //   // return *prev_it;
    //     // }
    //     // return (*it + supposed_pos) / 2;
    //     // // return *it;
    //   };

    //   // TODO: there's another solution: if the apperances of next kmer has
    //   for (const auto& e : out_edges) {
    //     // use pos as nearby is not accurate

    //     auto u = g.target(e);
    //     stk.emplace(u, pos + 1);

    //     break;
    //     // if (auto next_pos = is_nearby(g[u].appearances); next_pos != -1) {
    //     //   stk.emplace(u, next_pos);
    //     //   break;
    //     // } else {
    //     //   spdlog::debug("pos = {}, path.size() = {}, cnt = {}, {}({} - {})
    //     -> {}({} - {})", pos + 1, path.size(), g[e].count, g[v].kmer,
    //     //                 *g[v].appearances.begin(),
    //     *g[v].appearances.rbegin(),
    //     //                 g[g.target(e)].kmer,
    //     //                 *g[g.target(e)].appearances.begin(),
    //     //                 *g[g.target(e)].appearances.rbegin());
    //     // }
    //   }
    // }
    // return concat_vertices(path);
  }

  // auto print_graphviz(const Path& path) {
  //   auto p = fs::path("/mnt/ec/ness/yolkee/thesis/tests/01_10002.dot");
  //   std::ofstream fout(p);

  //   assert(fout.is_open());
  //   std::set<Vertex> vertexes;
  //   std::set<Edge> edges;
  //   for (const auto& [v, w] : path) {
  //     vertexes.insert(v);
  //   }

  //   std::function<void(Vertex)> dfs = [&](Vertex now) {
  //     vertexes.insert(now);
  //     auto sum = out_weight_sum(now);
  //     for (const auto& e : g.out_edges(now, false)) {
  //       auto u = g.target(e);
  //       if (g[e].count / sum > param.VALID_BRANCH_RATIO) {
  //         edges.insert(e);
  //         if (vertexes.contains(u)) {
  //           continue;
  //         }
  //         dfs(u);
  //       }
  //     }
  //   };

  //   for (const auto& v : vertexes) {
  //     auto sum = out_weight_sum(v);
  //     for (const auto& e : g.out_edges(v, false)) {
  //       auto u = g.target(e);
  //       if (g[e].count / sum > param.VALID_BRANCH_RATIO) {
  //         if (vertexes.contains(u)) {
  //           edges.insert(e);
  //         }
  //         dfs(u);
  //       }
  //     }
  //   }

  //   spdlog::debug("vertexes.size() = {}", vertexes.size());
  //   spdlog::debug("edges.size() = {}", edges.size());

  //   fout << "digraph G {\n";
  //   for (const auto& v : vertexes) {
  //     fout << "  " << v << " [label=\"" << g[v].kmer << "\"];\n";
  //   }
  //   for (const auto& e : edges) {
  //     fout << "  " << g.source(e) << " -> " << g.target(e) << " [label=\""
  //          << g[e].count << "\"];\n";
  //   }
  //   fout << "}\n";
  //   fout.close();

  //   /* transform dot to svg */
  //   auto svg = p;
  //   svg.replace_extension(".svg");
  //   // sfdp -x -Goverlap=scale -Tpng data.dot > data.png
  //   std::string cmd = fmt::format("sfdp -x -Goverlap=scale -Tsvg {} -o {}",
  //                                 p.string(), svg.string());
  //   spdlog::debug("start dot");
  //   std::system(cmd.c_str());
  //   spdlog::debug("end dot");
  // }

  auto print() {
    spdlog::debug("There are {} vertices.", g.vertices().size());
    spdlog::debug("There are {} kmers", kmer_to_vertex.size());
    auto start_addr = fmt::ptr(kmer_to_vertex.begin()->first.data());
    auto end_addr = fmt::ptr(kmer_to_vertex.rbegin()->first.data());
    spdlog::debug("Kmer start at {} and end at {}", start_addr, end_addr);
    spdlog::debug("Capacity = {}, size = {}", kmer_to_vertex.capacity(),
                  kmer_to_vertex.size());
    spdlog::debug("sz = {:.3f} MB", ((double) kmer_to_vertex.capacity() * sizeof(Kmer) / 1024 / 1024));

    // spdlog::debug("Size in memory = {}", sizeof(*this));
  }

  auto get_graph_size() const noexcept -> std::size_t {
    return std::ranges::ssize(g.vertices());
  }

  auto clear() {
    g.clear();
    kmer_to_vertex.clear();
    max_path.clear();
  }

  /* The underlying graph */
  Graph g;
  /* kmer size */
  std::size_t kmer_size;
  /* the array size that store the data, which equal to `kmer_size` /  */
  std::size_t used_bytes;
  std::size_t shift_at_last_unit;

  /* The maximum length that `get_read` can reach */
  std::size_t max_assemble_len;

  /* for debug, recording info of longest path */
  Path max_path;
  std::size_t max_path_size;

  /* each kmer and its vertex */
  KmerMap<Vertex> kmer_to_vertex;

  // debug
  std::map<Vertex, std::string> vertex_to_kmer;
  std::map<Edge, std::size_t> edge_to_weight;
};