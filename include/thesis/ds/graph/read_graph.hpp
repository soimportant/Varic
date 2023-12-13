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

#include <biovoltron/algo/assemble/graph/graph_wrapper.hpp>

namespace bio = biovoltron;

struct ReadGraph {
public:
  struct Param {

    /* the length that used to refind the source when not reaching sink */
    const std::size_t SEARCH_LEN_WHEN_NOT_GOAL = 100ul;

    /* The valid weight ratio that should consider as a branch */
    const double VALID_BRANCH_RATIO = 0.2;
  } param;

  struct VertexProperty {
    std::string kmer;
    std::set<std::size_t> appearances;
  };

  struct EdgeProperty {
    std::size_t count = 0;
  };

  using Graph = typename bio::GraphWrapper<VertexProperty, EdgeProperty>;
  using Vertex = typename Graph::Vertex;
  using Edge = typename Graph::Edge;
  using Path = std::vector<std::pair<Vertex, std::size_t>>;

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

  auto create_vertex(const std::string_view kmer) {
    // auto [iter, success] = unique_kmers.insert_or_assign(kmer, Vertex{});
    if (auto it = unique_kmers.find(kmer); it != unique_kmers.end()) {
      return it->second;
    }
    const auto v = g.create_vertex();
    g[v].kmer = kmer;
    unique_kmers[g[v].kmer] = v;
    return v;
  }

  auto create_edge(const Vertex &u, const Vertex &v) {
    const auto e = g.create_edge(u, v);
    g[e].count += 1;
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

  

public:
  ReadGraph(const std::size_t kmer_size, const std::size_t min_occ, const std::size_t max_assemble_len) {
    this->kmer_size = kmer_size;
    this->min_occ = min_occ;
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

  auto get_sources(const std::size_t start, const std::size_t end) {
    auto candidate_sources = get_vertex_inside_range(start, end);
    spdlog::debug("Find Sources({} - {}): candidate.size() = {}", start, end,
                  candidate_sources.size());
    for (auto min_weight = this->min_occ; min_weight >= 2; min_weight--) {
      spdlog::debug("searching for weight_threshold = {}...", min_weight);
      auto sources = std::vector<Vertex>{};
      std::ranges::copy_if(candidate_sources, std::back_inserter(sources),
                           [this, min_weight](const auto &v) {
                             return is_source(v, min_weight);
                           });
      if (sources.size() != 0) {
        std::ranges::sort(sources, [this](const auto &a, const auto &b) {
          const auto a_diff =
              *g[a].appearances.rbegin() - *g[a].appearances.begin();
          const auto b_diff =
              *g[b].appearances.rbegin() - *g[b].appearances.begin();
          return a_diff < b_diff;
          // return *g[a].appearances.rbegin() < *g[b].appearances.rbegin();
          // return std::tie(*g[a].appearances.begin(),
          // *g[a].appearances.rbegin()) <
          //        std::tie(*g[b].appearances.begin(),
          //        *g[b].appearances.rbegin());
          // return *g[a].appearances.begin() < *g[b].appearances.begin();
        });
        return sources;
      }
    }
    return std::vector<Vertex>{};
  }

  auto get_sinks(const std::size_t start, const std::size_t end) {
    auto candidate_sinks = get_vertex_inside_range(start, end);
    spdlog::debug("Find Sinks({} - {}): candidate.size() = {}", start, end,
                  candidate_sinks.size());

    for (auto min_weight = this->min_occ; min_weight >= 2; min_weight--) {
      spdlog::debug("searching for weight_threshold = {}...", min_weight);
      auto sinks = std::vector<Vertex>{};
      std::ranges::copy_if(
          candidate_sinks, std::back_inserter(sinks),
          [this, min_weight](const auto &v) { return is_sink(v, min_weight); });
      if (sinks.size() != 0) {
        std::ranges::sort(sinks, [this](const auto &a, const auto &b) {
          const auto a_diff =
              *g[a].appearances.rbegin() - *g[a].appearances.begin();
          const auto b_diff =
              *g[b].appearances.rbegin() - *g[b].appearances.begin();
          return a_diff < b_diff;
          // return *g[a].appearances.rbegin() < *g[b].appearances.rbegin();
          // return std::tie(*g[a].appearances.begin(),
          // *g[a].appearances.rbegin()) <
          //        std::tie(*g[b].appearances.begin(),
          //        *g[b].appearances.rbegin());
          // return *g[a].appearances.begin() < *g[b].appearances.begin();
        });
        return sinks;
      }
    }
    return std::vector<Vertex>{};
  }

  // auto path_finder(const Vertex& source, const Vertex& sink,
  //                  Path& path, std::map<Vertex, int>& vis,
  //                  std::set<Vertex>& fail_vertexes) {

  //   auto stk = std::stack<Vertex>{};
  //   auto inside_stk = std::set<Vertex>{};
  //   stk.push(source);
  //   inside_stk.insert(source);

  //   while (!stk.empty()) {
  //     auto v = stk.top();
  //     stk.pop();
  //     // inside_stk.erase(v);
  //     if (v == sink) {
  //       spdlog::debug("Reach sink {}", g[v].kmer);
  //       return true;
  //     }

  //     vis[v] += 1;
  //     auto out_edges = g.out_edges(v, false);
  //     std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b)
  //     {
  //       const auto a_vis_times = std::max(1, vis[g.target(a)]);
  //       const auto b_vis_times = std::max(1, vis[g.target(b)]);
  //       return static_cast<double>(g[a].count) / a_vis_times >
  //              static_cast<double>(g[b].count) / b_vis_times;
  //     });

  //     for (const auto& e : out_edges) {
  //       // if (g[e].count == 1) {
  //       //   continue;
  //       // }
  //       auto u = g.target(e);
  //       if (vis[u] >= g[e].count) {
  //         continue;
  //       }

  //       // if (inside_stk.contains(u)) {
  //       //   continue;
  //       // }
  //       // inside_stk.insert(u);

  //       stk.push(u);
  //       break;
  //     }
  //   }

  //   return false;

  //   // spdlog::debug("v = {}, ({} - {}), dep = {}", g[v].kmer,
  //   *g[v].appearances.begin(),
  //   //               *g[v].appearances.rbegin(), path.size());
  //   // // spdlog::debug("cnt = {}", cnt);
  //   // path.push_back(v);
  //   // if (v == sink) {
  //   //   spdlog::debug("Reach sink {}", g[v].kmer);
  //   //   cnt = 0;
  //   //   return true;
  //   // }
  //   // vis[v] += 1;
  //   // auto out_edges = g.out_edges(v, false);
  //   // std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto
  //   &b) {
  //   //   const auto a_vis_times = std::max(1, vis[g.target(a)]);
  //   //   const auto b_vis_times = std::max(1, vis[g.target(b)]);
  //   //   return static_cast<double>(g[a].count) / a_vis_times >
  //   //          static_cast<double>(g[b].count) / b_vis_times;
  //   // });
  //   // for (const auto& e : out_edges) {
  //   //   auto u = g.target(e);
  //   //   if (vis[u] >= g[e].count) {
  //   //     continue;
  //   //   }
  //   //   if (fail_vertexes.contains(u)) {
  //   //     continue;
  //   //   }
  //   //   if (path_finder(u, sink, path, vis, fail_vertexes)) {
  //   //     return true;
  //   //   }
  //   // }
  //   // cnt--;

  //   // fail_vertexes.insert(v);
  //   // path.pop_back();
  //   // return false;
  // }

  auto path_finder(const Vertex &now, const Vertex &goal, Path &path,
                   std::map<Vertex, std::size_t> &vis) {
    {
      max_path_size = std::max(max_path_size, path.size());
    }
    if (path.size() >= max_assemble_len) {
      spdlog::debug("Reach max_assemble_len {}", max_assemble_len);
      return false;
    }
    if (now == goal) {
      spdlog::debug("Reach sink {}", g[now].kmer);
      return true;
    }
    vis[now] += 1;
    auto out_edges = g.out_edges(now, false);
    auto edge_weights = std::vector<std::size_t>(out_edges.size());
    auto edge_weight_sum = 0.0;
    for (auto i = 0u; const auto &e : out_edges) {
      auto u = g.target(e);
      edge_weights[i++] = g[e].count / std::max(1ul, vis[u]);
      edge_weight_sum += edge_weights[i - 1];
    }
    auto order = std::vector<std::size_t>(out_edges.size());
    std::iota(order.begin(), order.end(), 0);
    std::ranges::sort(order, [&](const auto &a, const auto &b) {
      return edge_weights[a] > edge_weights[b];
    });
    for (const auto x : order) {
      auto [e, w] = std::tie(out_edges[x], edge_weights[x]);
      auto u = g.target(e);
      if (fail_vertexes.contains(u)) {
        continue;
      }
      if (w < param.VALID_BRANCH_RATIO * edge_weight_sum) {
        continue;
      }
      path.emplace_back(u, w);
      if (path_finder(u, goal, path, vis)) {
        return true;
      }
      path.pop_back();
    }
    fail_vertexes.insert(now);
    vis[now] -= 1;
    return false;
  }

  auto get_read(const Vertex &source, const Vertex &sink) {
    spdlog::debug("Try source = {}({} - {}) sink = {}({} - {})", g[source].kmer,
                  *g[source].appearances.begin(),
                  *g[source].appearances.rbegin(), g[sink].kmer,
                  *g[sink].appearances.begin(), *g[sink].appearances.rbegin());
    if (fail_vertexes.contains(source) || fail_vertexes.contains(sink)) {
      spdlog::debug("source or sink is failed", g[source].kmer);
      return std::string{};
    }
    if (source == sink) {
      spdlog::debug("source == sink");
      return std::string{};
    }

    /* store (Vertex, edge_weight) for further locally reassemble */
    auto path = std::vector<std::pair<Vertex, std::size_t>>{};
    auto vis = std::map<Vertex, std::size_t>{};
    max_path_size = 0;
    auto result = path_finder(source, sink, path, vis);
    {
      spdlog::debug("result = {}, max_path_size = {}", result, max_path_size);
    }
    // {
    //   for (auto [v, weight] : path) {
    //     spdlog::debug("v = {}, ({} - {}), weight = {}", g[v].kmer,
    //                   *g[v].appearances.begin(), *g[v].appearances.rbegin(),
    //                   weight);
    //   }
    // }
    auto corrected_read = g[source].kmer;
    for (const auto &[v, weight] : path) {
      corrected_read += g[v].kmer.back();
    }
    return corrected_read;

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

  void print() {
    auto vertices = g.vertices();
    spdlog::debug("vertices size = {}", std::ranges::ssize(vertices));
  }

  auto get_graph_size() const noexcept -> std::size_t {
    return std::ranges::ssize(g.vertices());
  }

  /* graph and its properties */
  Graph g;
  std::size_t kmer_size;
  std::size_t min_occ;

  /* The maximum length that `get_read` can reach*/
  std::size_t max_assemble_len;
  /* record the vertexes when failed on assembly */
  std::set<Vertex> fail_vertexes;

  /* for debug, recording  */
  std::size_t max_path_size;

  std::vector<std::size_t> order;

  // ? the usage of dup_kmers
  std::set<std::string> dup_kmers;

  /* each kmer and its vertex */
  std::map<std::string_view, Vertex, std::less<>> unique_kmers;
};