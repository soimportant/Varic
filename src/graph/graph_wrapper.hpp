#pragma once

#include <biovoltron/algo/align/inexact_match/smithwaterman.hpp>
#include <biovoltron/file_io/sam.hpp>
#include <biovoltron/utility/haplotype/haplotype.hpp>
#include <biovoltron/utility/read/quality_utils.hpp>

/* will remove it in future */
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <spdlog/spdlog.h>

namespace biovoltron {

template<class T>
auto boost_iterator_to_range(const std::pair<T, T>& p) {
  auto&& [beg, end] = p;
  return std::ranges::subrange(beg, end);
}

// TODO: has_cycle


/**
 * @brief Wrapper of Boost Grpah Libaray(BGL), the current implemetation is 
 * bidirecitonal graph
 * @todo add a template parameter for whether the graph is bidirecitonal or not.
 * @tparam VertexProperty 
 * @tparam EdgeProperty 
 */
template<class VertexProperty, class EdgeProperty> 
struct GraphWrapper {

public:
  using Graph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS,
    VertexProperty,
    EdgeProperty
  >;
  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = boost::graph_traits<Graph>::edge_descriptor;
  using Path = std::vector<Vertex>;

protected:
  using EdgePredicate = std::function<bool(Edge)>;
  EdgePredicate edge_filter = [](const auto&) { return true; };
  Graph g;

public:
  /**
   * @brief Create a vertex;
   * 
   * @return Vertex
   */
  auto create_vertex() {
    return boost::add_vertex(g);
  }

  /**
   * @brief Create an edge if not exists
   * 
   * @param u 
   * @param v 
   * @return Edge
   */
  auto create_edge(Vertex u, Vertex v) {
    const auto [e, success] = boost::add_edge(u, v, g);
    return e;
  }

  template<std::convertible_to<EdgePredicate> Predicate>
  auto set_edge_filter(Predicate&& filter) {
    edge_filter = filter;
  }


  /**
   * @brief find all paths between two vertices
   * 
   * @param from source Vertex
   * @param to destination Vertex
   * @return std::vector<Path>
   */
  auto find_paths(Vertex from, Vertex to, bool use_edge_filter = true) {
    auto paths = std::vector<Path>{};
    auto path = Path{};
    std::function<void(Vertex)> dfs = [&](Vertex u) {
      path.push_back(u);
      if (to == u) {
        paths.push_back(path);
      } else {
        for (auto e : out_edges(u, use_edge_filter)) {
          const auto v = target(e);
          if (std::ranges::find(path, v) == path.end()) {
            dfs(v);
          }
        }
      }
      path.pop_back();
    };
    dfs(from);
    return paths;
  }

  /* vertex properties */
  decltype(auto) operator[](const Vertex& v) {
    return (g[v]);
  }

  auto vertices() const {
    return boost_iterator_to_range(boost::vertices(g));
  }

  auto edges() const {
    return boost_iterator_to_range(boost::edges(g));
  }

  auto in_edges(const Vertex& v, const bool use_edge_filter = true) const {
    auto r = boost_iterator_to_range(boost::in_edges(v, g));
    std::vector<Edge> edges;
    if (use_edge_filter) {
      std::ranges::copy_if(r, std::back_inserter(edges), edge_filter);
    } else {
      std::ranges::copy(r, std::back_inserter(edges));
    }
    return edges;
  }

  auto out_edges(const Vertex& v, const bool use_edge_filter = true) const {
    auto r = boost_iterator_to_range(boost::out_edges(v, g));
    std::vector<Edge> edges;
    if (use_edge_filter) {
      std::ranges::copy_if(r, std::back_inserter(edges), edge_filter);
    } else {
      std::ranges::copy(r, std::back_inserter(edges));
    }
    return edges;
  }

  auto out_degree(const Vertex& v, const bool use_edge_filter = true) const {
    /**
     * ? Why not using std::ranges::size here?
     * std::views::filter_view doesn't meet the requirement to use
     * std::ranges::size(), which need the iterator type be random accessible,
     * so we use std::ranges::distance instead.
     */
    return std::ranges::distance(out_edges(v, use_edge_filter));
  }

  auto in_degree(const Vertex& v, const bool use_edge_filter = true) const {
    return std::ranges::distance(in_edges(v, use_edge_filter));
  }

  auto is_source(const Vertex& v, const bool use_edge_filter = true) const {
    return in_degree(v, use_edge_filter) == 0 && 
           out_degree(v, use_edge_filter) > 0;
  }

  auto is_sink(const Vertex& v, bool use_edge_filter = true) const {
    return out_degree(v, use_edge_filter) == 0 &&
           in_degree(v, use_edge_filter) > 0;
  }

  /* edge properties */
  decltype(auto) operator[](const Edge& e) {
    return (g[e]);
  }

  auto edge(const Vertex& u, const Vertex& v) const {
    /* return std::pair<Edge, bool> -> indicate whether edge is exist */
    return boost::edge(u, v, g).first;
  }

  auto source(const Edge& e) const {
    return boost::source(e, g);
  }

  auto target(const Edge& e) const {
    return boost::target(e, g);
  }

  /* graph properties */
  auto get_sources(const bool use_edge_filter = true) const {
    std::vector<Vertex> sources;
    std::ranges::copy_if(vertices(), std::back_inserter(sources), 
      [this, use_edge_filter](const Vertex &v) {
        return is_source(v, use_edge_filter);
      }
    );
    return sources;
  }

  auto get_sinks(const bool use_edge_filter = true) const {
    std::vector<Vertex> sinks;
    std::ranges::copy_if(vertices(), std::back_inserter(sinks),
      [this, use_edge_filter](const Vertex& v) {
        return is_sink(v, use_edge_filter);
      }
    );
    return sinks;
  }
};

// friend auto& operator<<(std::ostream& os, const GraphWrapper& wrapper) {
//   const auto& g = wrapper.g;
//   os << "digraph assembly_graphs {";
//   for (const auto e : boost::make_iterator_range(boost::edges(g))) {
//     os << boost::source(e, g) << " -> " << boost::target(e, g) << " ";
//     const auto count = g[e].count;
//     if (g[e].is_ref)
//       os << "[label=" << count << ",color=red];\n";
//     else if (count < PRUNE_FACTOR)
//       os << "[label=" << count << ",style=dotted,color=grey];\n";
//     else
//       os << "[label=" << count << "];\n";
//   }

//   for (const auto v : boost::make_iterator_range(boost::vertices(g))) {
//     os << v << " ";
//     const auto kmer = g[v].kmer;
//     if (boost::in_degree(v, g) == 0)
//       os << "[label=" << kmer << ",shape=box]\n";
//     else
//       os << "[label=" << kmer.back() << ",shape=box]\n";
//   }
//   return os << "}";
// }

} // namespace biovoltron