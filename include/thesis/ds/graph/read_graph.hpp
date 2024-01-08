#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <numeric>
#include <ranges>
#include <set>
#include <stack>
#include <string>
#include <string_view>
#include <vector>

#include <boost/graph/graphviz.hpp>

#include <biovoltron/algo/assemble/graph/graph_wrapper.hpp>

namespace bio = biovoltron;
namespace fs = std::filesystem;

struct ReadGraph {
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

  struct VertexProperty {
    std::string_view kmer;
    // std::string kmer;
    std::set<std::size_t> appearances;
  };

  struct EdgeProperty {
    std::size_t count = 0;
  };

  using Graph = typename bio::GraphWrapper<VertexProperty, EdgeProperty>;
  using Vertex = typename Graph::Vertex;
  using Edge = typename Graph::Edge;
  using Path = std::vector<std::pair<Vertex, std::size_t>>;

  auto create_vertex(const std::string_view kmer) {
    // auto [iter, success] = unique_kmers.insert_or_assign(kmer, Vertex{});
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) {
      return it->second;
    }
    const auto v = g.create_vertex();
    g[v].kmer = kmer;
    unique_kmers[kmer] = v;

    // debug
    vertex_to_kmer[v] = kmer;
    return v;
  }

  auto create_edge(const Vertex &u, const Vertex &v) {
    const auto e = g.create_edge(u, v);
    g[e].count += 1;

    // debug
    edge_to_weight[e] = g[e].count;
    return e;
  }

  auto get_vertex(const std::string_view kmer) {
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) {
      return it->second;
    }
    return create_vertex(kmer);
  }

  auto get_edge(const Vertex &u, const Vertex &v) {
    for (const auto e : g.out_edges(u, false)) {
      if (g.target(e) == v) {
        return e;
      }
    }
    return create_edge(u, v);
  }

  auto out_weight_sum(const Vertex &v) {
    auto out_edges = g.out_edges(v, false);
    return std::accumulate(out_edges.begin(), out_edges.end(), 0.0,
                           [this](const auto &a, const auto &b) {
                             return a + g[b].count;
                           });
  }

  auto in_weight_sum(const Vertex &v) {
    auto in_edges = g.in_edges(v, false);
    return std::accumulate(in_edges.begin(), in_edges.end(), 0.0,
                           [this](const auto &a, const auto &b) {
                             return a + g[b].count;
                           });
  }

  auto get_vertex_inside_range(const std::size_t start, const std::size_t end) {
    assert(start < end);
    auto vertices = std::vector<Vertex>{};
    for (const auto &[kmer, v] : unique_kmers) {
      auto it = g[v].appearances.lower_bound(start);
      if (it != g[v].appearances.end() && *it <= end) {
        vertices.push_back(v);
      }
    }
    return vertices;
  }

  auto is_source(const Vertex &v, const std::size_t weight_threshold) {
    auto out_degree =
        std::ranges::count_if(g.out_edges(v, false), [&, this](const auto &e) {
          return g[e].count >= weight_threshold;
        });
    auto in_degree =
        std::ranges::count_if(g.in_edges(v, false), [&, this](const auto &e) {
          return g[e].count >= weight_threshold;
        });
    return in_degree == 0 && out_degree > 0;
  }

  auto is_sink(const Vertex &v, const std::size_t weight_threshold) {
    auto out_degree =
        std::ranges::count_if(g.out_edges(v, false), [&, this](const auto &e) {
          return g[e].count >= weight_threshold;
        });
    auto in_degree =
        std::ranges::count_if(g.in_edges(v, false), [&, this](const auto &e) {
          return g[e].count >= weight_threshold;
        });
    return in_degree > 0 && out_degree == 0;
  }

  template<bool Reverse = false>
  auto concat_vertices(const Vertex& source, Path& path) {
    auto read = std::string(g[source].kmer);
    for (const auto& [v, weight] : path) {
      if constexpr (Reverse == false) {
        read += g[v].kmer.back();
      } else {
        read += g[v].kmer.front();
      }
    }
    if constexpr (Reverse) {
      std::ranges::reverse(read);
    }
    return read;
  }

  // TODO: extend source and sink
  template <bool Reverse = false> auto extend_path_at_end(Path &path) {}

  /* do local assemble between source and sink */
  template <bool Reverse = false>
  auto local_assemble(const Vertex &source, const Vertex &sink,
                      const std::size_t length_limit) {

    // spdlog::debug("Local assemble between {} - {}", g[source].kmer,
    //               g[sink].kmer);
    // TODO: considering following stuffs
    // 1. cycle
    // 2. the possibility of jumping to vertexes behind the sink
    //   - set a threshold of the path length
    auto paths = std::vector<Path>{};
    auto dfs = [&](auto &&self, const Vertex &now, Path &path,
                   std::set<Vertex> &vis) {
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
      for (const auto &e : edges) {
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
    /* the path doesn't contain source for preserving weight between source and its previous vertex */
    auto path = Path{};
    auto vis = std::set<Vertex>{};
    dfs(dfs, source, path, vis);
    // sort the paths by average weight
    std::ranges::sort(paths, [](const auto &a, const auto &b) {
      auto a_sum = std::accumulate(
          a.begin(), a.end(), 0.0,
          [](const auto &a, const auto &b) { return a + b.second; });
      auto b_sum = std::accumulate(
          b.begin(), b.end(), 0.0,
          [](const auto &a, const auto &b) { return a + b.second; });
      return a_sum / a.size() > b_sum / b.size();
    });
    // spdlog::debug("There are {} paths.", paths.size());
    for (auto &path : paths) {
      auto weight_sum = std::accumulate(
          path.begin(), path.end(), 0.0,
          [](const auto &a, const auto &b) { return a + b.second; });
      auto [mn, mx] = std::ranges::minmax_element(
          path.begin(), path.end(),
          [](const auto &a, const auto &b) { return a.second < b.second; });
      auto s = std::string{};
      for (auto [v, w] : path) {
        s.push_back(g[v].kmer.back());
      }
      // spdlog::debug("candidate({}, {}, {:.3f}): {}, ", mn->second, mx->second,
      //               weight_sum / path.size(), s);
    }
    assert(paths.size() > 0);
    return paths[0];
  }

  template<bool Reverse = false>
  auto valid_path(const Path& path) {
    auto sz = path.size();
    for (auto i = 1u; i < sz; i++) {
      auto [prev, prev_weight] = path[i - 1];
      auto [now, now_weight] = path[i];
      if constexpr (Reverse == false) {
        if (g[prev].kmer.substr(1) != g[now].kmer.substr(0, kmer_size - 1)) {
          return false;
        }
      } else {
        if (g[prev].kmer.substr(0, kmer_size - 1) !=
            g[now].kmer.substr(1)) {
          return false;
        }
      }
    }
    return true;
  }

  // check the path for the need of local reassemble
  template <bool Reverse = false>
  auto check_and_reassemble_path(const Path &path) {
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
        // spdlog::debug("Found branch point at {}, i = {}", g[path[i].first].kmer, i);
        auto prev_branch_pos = find_previous_breakpoint(path[i].second * param.LOCAL_ASSEMBLE_WEIGHT_RATIO);
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
              prev_branch_v, now_branch_v, old_path_len * param.LOCAL_ASSEMBLE_LENGTH_RATIO);
          // spdlog::debug("Found sink point at {}, len = {}", g[now_branch_v].kmer,
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
        // for (auto [v, w] : new_path | std::views::reverse | std::views::take(10) | std::views::reverse) {
        //   spdlog::debug("v = {}, weight = {}", g[v].kmer, w);
        // }
        // assert(valid_path<Reverse>(new_path));
      }
    }
    assert(valid_path<Reverse>(new_path));
    return new_path;

    // TODO: remove low weight path at the end
  }

  template<bool Reverse = false>
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
        // spdlog::debug("Reach sink {}, branch_cnt = {}", g[v].kmer, branch_cnt);
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
      for (auto i = 0u; const auto &e : edges) {
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
      std::ranges::sort(order, [&](const auto &a, const auto &b) {
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

  template <bool Reverse = false>
  auto path_finder(const Vertex &now, const Vertex &goal, Path &path,
                   std::map<Vertex, std::size_t> &vis,
                   std::set<Vertex> &fail_vertexes) {
    if (path.size() >= max_assemble_len) {
      // stuck in cycle, need increase kmer_size
      // spdlog::debug("Reach max_assemble_len {}", max_assemble_len);
      return false;
    }
    if (now == goal) {
      // spdlog::debug("Reach sink {}, path.size() = {}", g[now].kmer,
      //               path.size());
      return true;
    }
    vis[now] += 1;
    auto edges = std::vector<Edge>{};
    if constexpr (Reverse == false) {
      edges = g.out_edges(now, false);
    } else {
      edges = g.in_edges(now, false);
    }
    auto edge_weights = std::vector<std::size_t>(edges.size());
    auto edge_weight_sum = 0.0;
    for (auto i = 0u; const auto &e : edges) {
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
    std::ranges::sort(order, [&](const auto &a, const auto &b) {
      return edge_weights[a] > edge_weights[b];
    });
    for (const auto x : order) {
      auto [e, w] = std::tie(edges[x], edge_weights[x]);
      auto u = Vertex{};
      if constexpr (Reverse == false) {
        u = g.target(e);
      } else {
        u = g.source(e);
      }
      if (fail_vertexes.contains(u)) {
        continue;
      }
      if (w == 0 || w < param.VALID_BRANCH_RATIO * edge_weight_sum) {
        continue;
      }
      path.emplace_back(u, w);
      if (path_finder<Reverse>(u, goal, path, vis, fail_vertexes)) {
        return true;
      }
      path.pop_back();
    }
    vis[now] -= 1;
    fail_vertexes.insert(now);
    {
      if (path.size() > max_path_size) {
        max_path_size = path.size();
        /* record the path */
        max_path = path;
      }
    }
    return false;
  }

  auto reset_vertex_weight_on_path(Path &path) {
    // set the duplicated kmer weight to the minimum weight it appear
    auto min_weight = std::map<ReadGraph::Vertex, std::size_t>{};
    for (const auto &[vertex, weight] : path) {
      if (min_weight.contains(vertex)) {
        min_weight[vertex] = std::min(min_weight[vertex], weight);
      } else {
        min_weight[vertex] = weight;
      }
    }
    for (auto &[vertex, weight] : path) {
      weight = min_weight[vertex];
    }
  }

public:
  ReadGraph(const std::size_t kmer_size,
            const std::size_t max_assemble_len) {
    this->kmer_size = kmer_size;
    this->max_assemble_len = max_assemble_len;
  }

  /**
   * @brief add one sequence to the graph
   * @param sequence the sequence to be added
   * @param left_bound the left boundary of the sequence in the read
   */
  auto add_seq(const std::string_view sequence, const std::size_t left_bound) {
    if (sequence.size() < kmer_size) {
      return;
    }
    auto update_appearance = [&](const Vertex &v, const std::size_t pos) {
      g[v].appearances.insert(pos + left_bound);
    };
    auto prev_vertex = get_vertex(sequence.substr(0, kmer_size));
    update_appearance(prev_vertex, 0);
    for (auto i = 1u; i + kmer_size <= sequence.size(); i++) {
      bool found = false;
      for (const auto &edge : g.out_edges(prev_vertex, false)) {
        const auto u = g.target(edge);
        if (g[u].kmer.back() == sequence[i + kmer_size - 1]) {
          g[edge].count += 1;
          update_appearance(u, i);
          prev_vertex = u;
          found = true;
          break;
        }
      }
      if (!found) {
        const auto v = get_vertex(sequence.substr(i, kmer_size));
        update_appearance(v, i);
        create_edge(prev_vertex, v);
        prev_vertex = v;
      }
    }
  }

  auto get_sources(const std::size_t start, const std::size_t end, const std::size_t weight_threshold) {
    auto candidate_sources = get_vertex_inside_range(start, end);
    // spdlog::debug("Find Sources({} - {}) using weight = {}: candidate.size() = {}", start, end, weight_threshold,
    //               candidate_sources.size());
    auto sources = std::vector<Vertex>{};
    std::ranges::copy_if(candidate_sources, std::back_inserter(sources),
                         [this, weight_threshold](const auto &v) {
                           return is_source(v, weight_threshold);
                         });
    std::ranges::sort(sources, [this](const auto &a, const auto &b) {
      const auto [a_st, a_ed] =
          std::tie(*g[a].appearances.begin(), *g[a].appearances.rbegin());
      const auto [b_st, b_ed] =
          std::tie(*g[b].appearances.begin(), *g[b].appearances.rbegin());
      const auto a_diff = a_ed - a_st;
      const auto b_diff = b_ed - b_st;
      if (a_diff != b_diff) {
        return a_diff < b_diff;
      }
      if (a_ed != b_ed) {
        return a_ed < b_ed;
      }
      return b_st > a_st;
    });
    // std::ranges::sort(sources, [this](const auto &a, const auto &b) {
    //   const auto a_diff =
    //       *g[a].appearances.rbegin() - *g[a].appearances.begin();
    //   const auto b_diff =
    //       *g[b].appearances.rbegin() - *g[b].appearances.begin();
    //   return a_diff < b_diff;
    //   // return *g[a].appearances.rbegin() < *g[b].appearances.rbegin();
    //   // return std::tie(*g[a].appearances.begin(),
    //   // *g[a].appearances.rbegin()) <
    //   //        std::tie(*g[b].appearances.begin(),
    //   //        *g[b].appearances.rbegin());
    //   // return *g[a].appearances.begin() < *g[b].appearances.begin();
    // });
    return sources;
  }

  auto get_sinks(const std::size_t start, const std::size_t end, const std::size_t weight_threshold) {
    auto candidate_sinks = get_vertex_inside_range(start, end);
    // spdlog::debug("Find Sinks({} - {}) using weight = {}: candidate.size() = {}", start, end, weight_threshold,
    //               candidate_sinks.size());

    auto sinks = std::vector<Vertex>{};
    std::ranges::copy_if(candidate_sinks, std::back_inserter(sinks),
                         [this, weight_threshold](const auto &v) {
                           return is_sink(v, weight_threshold);
                         });
    std::ranges::sort(sinks, [this](const auto &a, const auto &b) {
      const auto [a_st, a_ed] =
          std::tie(*g[a].appearances.begin(), *g[a].appearances.rbegin());
      const auto [b_st, b_ed] =
          std::tie(*g[b].appearances.begin(), *g[b].appearances.rbegin());
      const auto a_diff = a_ed - a_st;
      const auto b_diff = b_ed - b_st;
      if (a_diff != b_diff) {
        return a_diff < b_diff;
      }
      if (a_st != b_st) {
        return a_st > b_st;
      }
      return a_ed < b_ed;
    });
    //   std::ranges::sort(sinks, [this](const auto &a, const auto &b) {
    //     const auto a_diff =
    //         *g[a].appearances.rbegin() - *g[a].appearances.begin();
    //     const auto b_diff =
    //         *g[b].appearances.rbegin() - *g[b].appearances.begin();
    //     return a_diff < b_diff;
    //     // return *g[a].appearances.rbegin() < *g[b].appearances.rbegin();
    //     // return std::tie(*g[a].appearances.begin(),
    //     // *g[a].appearances.rbegin()) <
    //     //        std::tie(*g[b].appearances.begin(),
    //     //        *g[b].appearances.rbegin());
    //     // return *g[a].appearances.begin() < *g[b].appearances.begin();
    //   });
    //   return sinks;
    return sinks;
  }

  /**
   * @brief Get the read object
   * 
   * @param source 
   * @param sink ideal sink
   * @return auto 
   */
  template<bool Reverse = false>
  auto find_read(const Vertex &source, const Vertex &sink) {
    assert(source != sink);
    // spdlog::debug("Find the path in {} direction", Reverse ? "reverse" : "forward");
    // spdlog::debug("Try source = {}({} - {}) sink = {}({} - {})", g[source].kmer,
    //               *g[source].appearances.begin(),
    //               *g[source].appearances.rbegin(), g[sink].kmer,
    //               *g[sink].appearances.begin(), *g[sink].appearances.rbegin());

    /* reset max_path */
    max_path_size = 0;
    max_path.clear();

    /* store (Vertex, edge_weight) for further locally reassemble */
    auto path = Path{};
    
    // auto vis = std::map<Vertex, std::size_t>{};
    // auto fail_vertexes = std::set<Vertex>{};
    // auto result = path_finder<Reverse>(source, sink, path, vis, fail_vertexes);
    // if (result == false) {
    //   if (max_path_size == max_assemble_len) {
    //     // spdlog::debug("max_path length reach threshold {}", max_assemble_len);
    //     return concat_vertices<Reverse>(source, max_path);
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

  auto print_graphviz(const Path& path) {

    auto p = fs::path("/mnt/ec/ness/yolkee/thesis/tests/01_10002.dot");
    std::ofstream fout(p);

    assert(fout.is_open());
    std::set<Vertex> vertexes;
    std::set<Edge> edges;
    for (const auto& [v, w] : path) {
      vertexes.insert(v);
    }

    std::function<void(Vertex)> dfs = [&](Vertex now) {
      vertexes.insert(now);
      auto sum = out_weight_sum(now);
      for (const auto& e : g.out_edges(now, false)) {
        auto u = g.target(e);
        if (g[e].count / sum > param.VALID_BRANCH_RATIO) {
          edges.insert(e);
          if (vertexes.contains(u)) {
            continue;
          }
          dfs(u);
        }
      }
    };

    for (const auto& v : vertexes) {
      auto sum = out_weight_sum(v);
      for (const auto& e : g.out_edges(v, false)) {
        auto u = g.target(e);
        if (g[e].count / sum > param.VALID_BRANCH_RATIO) {
          if (vertexes.contains(u)) {
            edges.insert(e);
          }
          dfs(u);
        }
      }
    }

    spdlog::debug("vertexes.size() = {}", vertexes.size());
    spdlog::debug("edges.size() = {}", edges.size());

    fout << "digraph G {\n";
    for (const auto &v : vertexes) {
      fout << "  " << v << " [label=\"" << g[v].kmer << "\"];\n";
    }
    for (const auto &e : edges) {
      fout << "  " << g.source(e) << " -> " << g.target(e)
                << " [label=\"" << g[e].count << "\"];\n";
    }
    fout << "}\n";
    fout.close();

    /* transform dot to svg */
    auto svg = p;
    svg.replace_extension(".svg");
    // sfdp -x -Goverlap=scale -Tpng data.dot > data.png
    std::string cmd =
        fmt::format("sfdp -x -Goverlap=scale -Tsvg {} -o {}", p.string(), svg.string());
    spdlog::debug("start dot");
    std::system(cmd.c_str());
    spdlog::debug("end dot");
  }

  auto get_graph_size() const noexcept -> std::size_t {
    return std::ranges::ssize(g.vertices());
  }

  auto clear() {
    g.clear();
    unique_kmers.clear();
    max_path.clear();
  }

  /* graph and its properties */
  Graph g;
  std::size_t kmer_size;

  /* The maximum length that `get_read` can reach*/
  std::size_t max_assemble_len;

  /* for debug, recording info of longest path */
  Path max_path;
  std::size_t max_path_size;


  /* each kmer and its vertex */
  std::map<std::string_view, Vertex, std::less<>> unique_kmers;

  // debug
  std::map<Vertex, std::string_view> vertex_to_kmer;
  std::map<Edge, std::size_t> edge_to_weight;
};