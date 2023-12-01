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
    const std::size_t MAX_OCCURANCE_AT_END = 4ul;
    std::size_t NEARBY_THRESHOLD;
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
  using Path = std::vector<Vertex>;

  /* graph and its properties */
  Graph g;
  std::size_t kmer_size;
  std::size_t minimum_occurance;

  // ? the usage of dup_kmers
  std::set<std::string> dup_kmers;

  /* each kmer and its vertex */
  std::map<std::string_view, Vertex, std::less<>> unique_kmers;

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

  auto concat_vertices(const Path &path) {
    if (path.empty()) {
      return std::string{};
    }
    auto u = path[0];
    auto seq = g[u].kmer;
    for (auto i = 1u; i < path.size(); i++) {
      auto v = path[i];
      seq += g[v].kmer.back();
    }
    return seq;
  }

public:
  ReadGraph(const std::size_t kmer_size, const std::size_t minimum_occurance)
      : kmer_size(kmer_size), minimum_occurance(minimum_occurance) {
    g.set_edge_filter([this](const Edge &e) {
      return g[e].count >= this->minimum_occurance;
    });
  }

  auto set_min_occ(const std::size_t min_occ) {
    minimum_occurance = min_occ;
    g.set_edge_filter([this](const Edge &e) {
      return g[e].count >= this->minimum_occurance;
    });
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

  auto get_sources(const std::size_t head_boundary) {
    auto possible_sources = std::vector<Vertex>{};
    for (const auto &[kmer, v] : unique_kmers) {
      auto it = g[v].appearances.begin();
      if (it != g[v].appearances.end() && *it <= head_boundary - kmer_size) {
        possible_sources.push_back(v);
      }
    }
    const auto threshold = std::min(
      param.MAX_OCCURANCE_AT_END,
      static_cast<std::size_t>(std::ceil(possible_sources.size() / 100.0))
    );
    spdlog::debug("head_len = {}, possible_sources.size() = {}", head_boundary, possible_sources.size());
    auto sources = std::vector<Vertex>{};
    for (const auto& source : possible_sources) {
      auto out_degree = std::ranges::count_if(
        g.out_edges(source, false), 
        [&, this](const auto& e) {
          return g[e].count >= threshold;
        });
      auto in_degree = std::ranges::count_if(
        g.in_edges(source, false), [&, this](const auto& e) {
          return g[e].count >= threshold;
        });
      if (in_degree == 0 && out_degree > 0) {
        sources.push_back(source);
      }
    }
    std::ranges::sort(sources, [this](const auto &a, const auto &b) {
      const auto a_diff = *g[a].appearances.rbegin() - *g[a].appearances.begin();
      const auto b_diff = *g[b].appearances.rbegin() - *g[b].appearances.begin();
      return a_diff < b_diff;
      // return *g[a].appearances.rbegin() < *g[b].appearances.rbegin();
      // return std::tie(*g[a].appearances.begin(), *g[a].appearances.rbegin()) <
      //        std::tie(*g[b].appearances.begin(), *g[b].appearances.rbegin());
      // return *g[a].appearances.begin() < *g[b].appearances.begin();
    });
    return sources;
  }

  auto get_sinks(const std::size_t tail_boundary) {
    auto possible_sinks = std::vector<Vertex>{};
    for (const auto& [kmer, v] : unique_kmers) {
      auto it = g[v].appearances.rbegin();
      if (it != g[v].appearances.rend() && *it >= tail_boundary) {
        possible_sinks.push_back(v);
      }
    }
    const auto threshold = std::min(
        param.MAX_OCCURANCE_AT_END,
        static_cast<std::size_t>(std::ceil(possible_sinks.size() / 100.0)));
    spdlog::debug("tail_len = {}, possible_sinks.size() = {}", tail_boundary, possible_sinks.size());
    auto sinks = std::vector<Vertex>{};
    for (const auto& sink : possible_sinks) {
      auto out_degree = std::ranges::count_if(
          g.out_edges(sink, false),
          [&, this](const auto &e) { return g[e].count >= threshold; });
      auto in_degree = std::ranges::count_if(
          g.in_edges(sink, false),
          [&, this](const auto &e) { return g[e].count >= threshold; });
      if (out_degree == 0 && in_degree > 0) {
        sinks.push_back(sink);
      }
    }
    std::ranges::sort(sinks, [this](const auto &a, const auto &b) {
      const auto a_diff =
          *g[a].appearances.rbegin() - *g[a].appearances.begin();
      const auto b_diff =
          *g[b].appearances.rbegin() - *g[b].appearances.begin();
      return a_diff < b_diff;
      // return *g[a].appearances.rbegin() > *g[b].appearances.rbegin();
      // return std::tie(*g[a].appearances.begin(), *g[a].appearances.rbegin()) >
      //        std::tie(*g[b].appearances.begin(), *g[b].appearances.rbegin());
    });
    return sinks;
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
  //     std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b) {
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



  //   // spdlog::debug("v = {}, ({} - {}), dep = {}", g[v].kmer, *g[v].appearances.begin(),
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
  //   // std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b) {
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

  auto get_read(const Vertex &source, const Vertex &sink, const std::size_t max_len) {
    spdlog::debug("Try source = {}({} - {}) sink = {}({} - {})", g[source].kmer,
                  *g[source].appearances.begin(),
                  *g[source].appearances.rbegin(), g[sink].kmer,
                  *g[sink].appearances.begin(), *g[sink].appearances.rbegin());
    auto path = Path{};
    auto vis = std::map<Vertex, int>{};
    auto fail_vertexes = std::set<Vertex>{};

    auto stk = std::stack<Vertex>{};
    stk.push(source);

    while (!stk.empty()) {
      auto v = stk.top();
      stk.pop();
      // spdlog::debug("v = {}, ({} - {}), dep = {}", g[v].kmer, *g[v].appearances.begin(),
      //               *g[v].appearances.rbegin(), path.size());
      // spdlog::debug("cnt = {}", cnt);
      path.push_back(v);
      /* if path.size() >= max_len -> stuck in cycle */
      if (v == sink || path.size() >= max_len) {
        // spdlog::debug("Reach sink {}", g[v].kmer);
        break;
      }
      vis[v] += 1;
      auto out_edges = g.out_edges(v, false);
      std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b) {
        const auto a_vis_times = std::max(1, vis[g.target(a)]);
        const auto b_vis_times = std::max(1, vis[g.target(b)]);
        return static_cast<double>(g[a].count) / a_vis_times >
               static_cast<double>(g[b].count) / b_vis_times;
      });
      for (const auto& e : out_edges) {
        auto u = g.target(e);
        if (vis[u] >= g[e].count) {
          continue;
        }
        stk.push(u);
        break;
      }
    }

    return concat_vertices(path);



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
    //   std::ranges::sort(out_edges, [this, &vis](const auto &a, const auto &b) {
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
    //     //   spdlog::debug("pos = {}, path.size() = {}, cnt = {}, {}({} - {}) -> {}({} - {})", pos + 1, path.size(), g[e].count, g[v].kmer,
    //     //                 *g[v].appearances.begin(), *g[v].appearances.rbegin(),
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
};