//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

/*
  This file implements the function

  template <class VertexAndEdgeListGraph, class DistanceMatrix,
            class P, class T, class R>
  bool
  johnson_all_pairs_shortest_paths
    (VertexAndEdgeListGraph& g, 
     DistanceMatrix& D,
     const bgl_named_params<P, T, R>& params)
 */

#ifndef BOOST_GRAPH_JOHNSON_HPP
#define BOOST_GRAPH_JOHNSON_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/ct_if.hpp>
#include <boost/type_traits/same_traits.hpp>

namespace boost {

  /*template <class VertexAndEdgeListGraph, class DistanceMatrix,
            class VertexID, class Weight, typename BinaryPredicate, 
            typename BinaryFunction, typename Infinity, class DistanceZero>
  bool
  johnson_all_pairs_shortest_paths(VertexAndEdgeListGraph& g1, 
               DistanceMatrix& D,
               VertexID id1, Weight w1, const BinaryPredicate& compare, 
               const BinaryFunction& combine, const Infinity& inf,
               DistanceZero zero)
  {
    typedef graph_traits<VertexAndEdgeListGraph> Traits1;
    typedef typename property_traits<Weight>::value_type DT;
    function_requires< BasicMatrixConcept<DistanceMatrix,
      typename Traits1::vertices_size_type, DT> >();

    typedef typename Traits1::directed_category DirCat;
    bool is_undirected = is_same<DirCat, undirected_tag>::value;

    typedef adjacency_list<vecS, vecS, directedS, 
      property< vertex_distance_t, DT>,
      property< edge_weight_t, DT, 
      property< edge_weight2_t, DT > > > Graph2;
    typedef graph_traits<Graph2> Traits2;

    Graph2 g2(num_vertices(g1) + 1);
    typename property_map<Graph2, edge_weight_t>::type 
      w = get(edge_weight, g2);
    typename property_map<Graph2, edge_weight2_t>::type 
      w_hat = get(edge_weight2, g2);
    typename property_map<Graph2, vertex_distance_t>::type 
      d = get(vertex_distance, g2);
    typedef typename property_map<Graph2, vertex_index_t>::type VertexID2;
    VertexID2 id2 = get(vertex_index, g2);

    // Construct g2 where V[g2] = V[g1] U {s}
    //   and  E[g2] = E[g1] U {(s,v)| v in V[g1]}
    std::vector<typename Traits1::vertex_descriptor> 
      verts1(num_vertices(g1) + 1);
    typename Traits2::vertex_descriptor s = *vertices(g2).first;
    {
      typename Traits1::vertex_iterator v, v_end;
      int i = 1;
      for (tie(v, v_end) = vertices(g1); v != v_end; ++v, ++i) {
        typename Traits2::edge_descriptor e; bool z;
        tie(e, z) = add_edge(s, get(id1, *v) + 1, g2);
        put(w, e, zero);
        verts1[i] = *v;
      }
      typename Traits1::edge_iterator e, e_end;
      for (tie(e, e_end) = edges(g1); e != e_end; ++e) {
        typename Traits2::edge_descriptor e2; bool z;
        tie(e2, z) = add_edge(get(id1, source(*e, g1)) + 1, 
                              get(id1, target(*e, g1)) + 1, g2);
        put(w, e2, get(w1, *e));
        if (is_undirected) {
          tie(e2, z) = add_edge(get(id1, target(*e, g1)) + 1, 
                                get(id1, source(*e, g1)) + 1, g2);
          put(w, e2, get(w1, *e));
        }
      }
    }
    typename Traits2::vertex_iterator v, v_end, u, u_end;
    typename Traits2::edge_iterator e, e_end;
    std::vector<DT> h_vec(num_vertices(g2));
    typedef typename std::vector<DT>::iterator iter_t;
    iterator_property_map<iter_t,VertexID2,DT,DT&> h(h_vec.begin(), id2);

    for (tie(v, v_end) = vertices(g2); v != v_end; ++v)
      d[*v] = inf;

    put(d, s, zero);
    // Using the non-named parameter versions of bellman_ford and
    // dijkstra for portability reasons.
    dummy_property_map pred; bellman_visitor<> bvis;
    if (bellman_ford_shortest_paths
        (g2, num_vertices(g2), w, pred, d, combine, compare, bvis)) {
      for (tie(v, v_end) = vertices(g2); v != v_end; ++v)
        put(h, *v, get(d, *v));
      // Reweight the edges to remove negatives
      for (tie(e, e_end) = edges(g2); e != e_end; ++e) {
        typename Traits2::vertex_descriptor a = source(*e, g2),
          b = target(*e, g2);
        put(w_hat, *e, get(w, *e) + get(h, a) - get(h, b));
      }
      for (tie(u, u_end) = vertices(g2); u != u_end; ++u) {
        dijkstra_visitor<> dvis;
        dijkstra_shortest_paths
          (g2, *u, pred, d, w_hat, id2, compare, combine, inf, zero,dvis);
        for (tie(v, v_end) = vertices(g2); v != v_end; ++v) {
          if (*u != s && *v != s) {
            typename Traits1::vertex_descriptor u1, v1;
            u1 = verts1[id2[*u]]; v1 = verts1[id2[*v]];
            D[id2[*u]-1][id2[*v]-1] = get(d, *v) + get(h, *v) - get(h, *u);
          }
        }
      }
      return true;
    } else
      return false;
  }*/

  namespace detail {
    template <class VertexAndEdgeListGraph, class VertexIDMap, 
      class WeightMap, typename DistanceZero>
    struct single_source_vertex_graph {
      VertexAndEdgeListGraph& _g;
      VertexIDMap _id;
      WeightMap _w;
      DistanceZero _zero;
      single_source_vertex_graph(VertexAndEdgeListGraph& g, 
        VertexIDMap id, WeightMap w, DistanceZero zero) 
        : _g(g), _id(id), _w(w), _zero(zero) {}
    };
    template <class Vertex, class Weight>
    struct ssv_graph_edge
        : std::pair<Vertex,Vertex> {
      Weight _w;
      ssv_graph_edge() {}
      ssv_graph_edge(Vertex s, Vertex t, Weight w) : 
        std::pair<Vertex,Vertex>(s,t), _w(w) {}
      //ssv_graph_edge(const ssv_graph_edge& e): _w(e._w) {}
    };
    template <class VertexAndEdgeListGraph, class EdgeIter, 
      class VertexIDMap, class WeightMap, class Edge, typename DistanceZero>
    class ssv_edge_iter :
      public iterator_facade<
        ssv_edge_iter<VertexAndEdgeListGraph,EdgeIter,
          VertexIDMap,WeightMap,Edge,DistanceZero>,
        Edge,
        forward_traversal_tag,
        Edge> 
    {
      typedef typename property_traits<VertexIDMap>::value_type src_iter;
      bool _iter_source;
      EdgeIter _ei, _eiend;
      src_iter _si,_siend;
      VertexAndEdgeListGraph *_g;
      VertexIDMap _id;
      WeightMap _w;
      DistanceZero _zero;
    public:
      ssv_edge_iter() {}
      ssv_edge_iter(bool iter_source, src_iter si, src_iter siend,
        EdgeIter ei, EdgeIter eiend, VertexAndEdgeListGraph& g, 
        VertexIDMap id, WeightMap w, DistanceZero zero)
        : _iter_source(iter_source), _si(si), _siend(siend), 
        _ei(ei), _eiend(ei), _g(&g), _id(id), _w(w), _zero(zero) {}
      ssv_edge_iter(const ssv_edge_iter& sei)
          : _iter_source(sei._iter_source), _si(sei._si), _siend(sei._siend),
          _ei(sei._ei), _eiend(sei._eiend), _g(sei._g), _id(sei._id),
          _w(sei._w), _zero(sei._zero) {}

    private:
      friend class boost::iterator_core_access;
      inline void increment() {
          if (_iter_source) { 
              ++_si; 
              if ( _si == _siend) { _iter_source = false; } 
          } else if (_ei != _eiend) { 
              ++_ei; 
          }
      }
      inline Edge dereference() const { 
          if (_iter_source) { return Edge(0,_si,_zero); }
          else { return Edge(_id[source(*_ei,*_g)]+1,_id[target(*_ei,*_g)]+1,_w[*_ei]); }
      }
      inline bool equal(const ssv_edge_iter& other) const {
          return (_iter_source == other._iter_source && 
                 _si == other._si && _ei == other._ei);
      }
    };
    template <class VertexAndEdgeListGraph, class VertexIDMap, 
      class WeightMap, typename DistanceZero>
    struct ssv_graph_help {
      typedef VertexAndEdgeListGraph native_graph;
      typedef graph_traits<VertexAndEdgeListGraph> native_traits;
      typedef typename native_traits::edge_iterator native_edge_iterator;
      typedef detail::ssv_graph_edge<
        typename property_traits<VertexIDMap>::value_type,
        typename property_traits<WeightMap>::value_type > edge;
      typedef typename ssv_edge_iter<native_graph, native_edge_iterator, 
          VertexIDMap, WeightMap, edge, DistanceZero> edge_iterator;
      typedef std::pair<edge_iterator, edge_iterator> edges_ret_type;
      typedef typename native_traits::edges_size_type num_edges_ret_type;
    };
    template <class Vertex, class Weight>
    class ssv_edge_weight_map {};
  }

  template <class VertexAndEdgeListGraph, class VertexIDMap, 
    class WeightMap, typename DistanceZero>
  struct graph_traits< detail::single_source_vertex_graph<
      VertexAndEdgeListGraph, VertexIDMap, WeightMap, DistanceZero> > {
    typedef edge_list_graph_tag traversal_category;
    typedef directed_tag directed_category;
    typedef typename graph_traits<VertexAndEdgeListGraph>::edge_parallel_category
        edge_parallel_category;
    typedef typename graph_traits<VertexAndEdgeListGraph>::edges_size_type
        edges_size_type;
    typedef typename property_traits<VertexIDMap>::value_type vertex_descriptor;
    typedef detail::ssv_graph_edge<vertex_descriptor, 
        typename property_traits<WeightMap>::value_type> edge_descriptor;
    typedef detail::ssv_graph_help<VertexAndEdgeListGraph,
      VertexIDMap, WeightMap, DistanceZero> helper;
    typedef typename helper::edge_iterator edge_iterator;
    static vertex_descriptor null_vertex() { 
        return std::numeric_limits<vertex_descriptor>::max(); }
  };

  /*template <class Vertex, class Weight, class Graph>
  Vertex source(detail::ssv_graph_edge<Vertex,Weight> e, Graph&) 
  { return e._s; }
  template <class Vertex, class Weight, class Graph>
  Vertex target(detail::ssv_graph_edge<Vertex,Weight> e, Graph&) 
  { return e._t; }*/
  template <class VertexAndEdgeListGraph, class VertexIDMap, 
    class WeightMap, typename DistanceZero>
  typename detail::ssv_graph_help<
    VertexAndEdgeListGraph, 
    VertexIDMap, WeightMap, DistanceZero>::edges_ret_type
  edges(const detail::single_source_vertex_graph<
          VertexAndEdgeListGraph, VertexIDMap, WeightMap, DistanceZero>& g) {
    typedef detail::single_source_vertex_graph<
          VertexAndEdgeListGraph, VertexIDMap, WeightMap, DistanceZero> graph;
    typedef typename graph_traits<VertexAndEdgeListGraph>::edge_iterator native_ei;
    native_ei ei,eiend;
    tie(ei,eiend) = edges(g._g);
    return std::make_pair(
        graph_traits<graph>::edge_iterator(true, 1, num_vertices(g._g)+1,
            ei, eiend, g._g, g._id, g._w, g._zero),
        graph_traits<graph>::edge_iterator(false, num_vertices(g._g)+1,
            num_vertices(g._g)+1, eiend, eiend, g._g, g._id, g._w, g._zero));
  }
  template <class VertexAndEdgeListGraph, class VertexIDMap, 
    class WeightMap, typename DistanceZero>
  typename detail::ssv_graph_help<
    VertexAndEdgeListGraph, VertexIDMap, 
      WeightMap, DistanceZero>::num_edges_ret_type
  num_edges(const detail::single_source_vertex_graph<
              VertexAndEdgeListGraph, VertexIDMap, WeightMap, DistanceZero>& g) {
    return num_edges(g._g) + num_vertices(g._g);
  }

  template <class Vertex, class Weight>
  Weight get(const detail::ssv_edge_weight_map<Vertex,Weight>&, 
    detail::ssv_graph_edge<Vertex,Weight>& e)
  { return e._w; }

  template <class Vertex, class Weight>
  struct property_traits<detail::ssv_edge_weight_map<Vertex,Weight> > {
    typedef readable_property_map_tag category;
    typedef Weight value_type;
    typedef detail::ssv_graph_edge<Vertex,Weight> key_type;
    typedef value_type reference;
  };

  template <class VertexAndEdgeListGraph, class DistanceMatrix,
            class VertexID, class Weight, typename BinaryPredicate, 
            typename BinaryFunction, typename Infinity, class DistanceZero>
  bool
  johnson_all_pairs_shortest_paths(VertexAndEdgeListGraph& g1, 
               DistanceMatrix& D, 
               VertexID id1, Weight w1, const BinaryPredicate& compare, 
               const BinaryFunction& combine, const Infinity& inf,
               DistanceZero zero)
  {
    typedef detail::single_source_vertex_graph<VertexAndEdgeListGraph,
        VertexID,Weight,DistanceZero> graph_ssv;

    typedef typename property_traits<Weight>::value_type DT;

    graph_ssv gssv(g1, id1, w1, zero);

    std::vector<DT> d_vec(num_vertices(g1)+1, inf);
    std::vector<DT> h_vec(num_vertices(g1)+1);
    typedef typename std::vector<DT>::iterator iter_t;
    iterator_property_map<iter_t,identity_property_map,DT,DT&> 
        d(d_vec.begin(), identity_property_map());
    iterator_property_map<iter_t,identity_property_map,DT,DT&> 
        h(h_vec.begin(), identity_property_map());

    detail::ssv_edge_weight_map<graph_traits<graph_ssv>::vertex_descriptor,
        property_traits<Weight>::value_type> w2;

    dummy_property_map pred; bellman_visitor<> bvis;
    if (bellman_ford_shortest_paths
      (gssv, num_vertices(g1)+1, w2, pred, d, combine, compare, bvis))
    {
        return true;
    } else
      return false;
    /*typedef graph_traits<VertexAndEdgeListGraph> Traits1;
    typedef typename property_traits<Weight>::value_type DT;
    function_requires< BasicMatrixConcept<DistanceMatrix,
      typename Traits1::vertices_size_type, DT> >();

    typedef typename Traits1::directed_category DirCat;
    bool is_undirected = is_same<DirCat, undirected_tag>::value;

    typedef adjacency_list<vecS, vecS, directedS, 
      property< vertex_distance_t, DT>,
      property< edge_weight_t, DT, 
      property< edge_weight2_t, DT > > > Graph2;
    typedef graph_traits<Graph2> Traits2;

    Graph2 g2(num_vertices(g1) + 1);
    typename property_map<Graph2, edge_weight_t>::type 
      w = get(edge_weight, g2);
    typename property_map<Graph2, edge_weight2_t>::type 
      w_hat = get(edge_weight2, g2);
    typename property_map<Graph2, vertex_distance_t>::type 
      d = get(vertex_distance, g2);
    typedef typename property_map<Graph2, vertex_index_t>::type VertexID2;
    VertexID2 id2 = get(vertex_index, g2);

    // Construct g2 where V[g2] = V[g1] U {s}
    //   and  E[g2] = E[g1] U {(s,v)| v in V[g1]}
    std::vector<typename Traits1::vertex_descriptor> 
      verts1(num_vertices(g1) + 1);
    typename Traits2::vertex_descriptor s = *vertices(g2).first;
    {
      typename Traits1::vertex_iterator v, v_end;
      int i = 1;
      for (tie(v, v_end) = vertices(g1); v != v_end; ++v, ++i) {
        typename Traits2::edge_descriptor e; bool z;
        tie(e, z) = add_edge(s, get(id1, *v) + 1, g2);
        put(w, e, zero);
        verts1[i] = *v;
      }
      typename Traits1::edge_iterator e, e_end;
      for (tie(e, e_end) = edges(g1); e != e_end; ++e) {
        typename Traits2::edge_descriptor e2; bool z;
        tie(e2, z) = add_edge(get(id1, source(*e, g1)) + 1, 
                              get(id1, target(*e, g1)) + 1, g2);
        put(w, e2, get(w1, *e));
        if (is_undirected) {
          tie(e2, z) = add_edge(get(id1, target(*e, g1)) + 1, 
                                get(id1, source(*e, g1)) + 1, g2);
          put(w, e2, get(w1, *e));
        }
      }
    }
    typename Traits2::vertex_iterator v, v_end, u, u_end;
    typename Traits2::edge_iterator e, e_end;
    std::vector<DT> h_vec(num_vertices(g2));
    typedef typename std::vector<DT>::iterator iter_t;
    iterator_property_map<iter_t,VertexID2,DT,DT&> h(h_vec.begin(), id2);

    for (tie(v, v_end) = vertices(g2); v != v_end; ++v)
      d[*v] = inf;

    put(d, s, zero);
    // Using the non-named parameter versions of bellman_ford and
    // dijkstra for portability reasons.
    dummy_property_map pred; bellman_visitor<> bvis;
    if (bellman_ford_shortest_paths
        (g2, num_vertices(g2), w, pred, d, combine, compare, bvis)) {
      for (tie(v, v_end) = vertices(g2); v != v_end; ++v)
        put(h, *v, get(d, *v));
      // Reweight the edges to remove negatives
      for (tie(e, e_end) = edges(g2); e != e_end; ++e) {
        typename Traits2::vertex_descriptor a = source(*e, g2),
          b = target(*e, g2);
        put(w_hat, *e, get(w, *e) + get(h, a) - get(h, b));
      }
      for (tie(u, u_end) = vertices(g2); u != u_end; ++u) {
        dijkstra_visitor<> dvis;
        dijkstra_shortest_paths
          (g2, *u, pred, d, w_hat, id2, compare, combine, inf, zero,dvis);
        for (tie(v, v_end) = vertices(g2); v != v_end; ++v) {
          if (*u != s && *v != s) {
            typename Traits1::vertex_descriptor u1, v1;
            u1 = verts1[id2[*u]]; v1 = verts1[id2[*v]];
            D[id2[*u]-1][id2[*v]-1] = get(d, *v) + get(h, *v) - get(h, *u);
          }
        }
      }
      return true;
    } else
      return false;*/
  }

  template <class VertexAndEdgeListGraph, class DistanceMatrix,
            class VertexID, class Weight, class DistanceZero>
  bool
  johnson_all_pairs_shortest_paths(VertexAndEdgeListGraph& g1, 
               DistanceMatrix& D,
               VertexID id1, Weight w1, DistanceZero zero)
  {
    typedef typename property_traits<Weight>::value_type WT;
    return johnson_all_pairs_shortest_paths(g1, D, id1, w1, 
                                            std::less<WT>(),
                                            closed_plus<WT>(),
                                            (std::numeric_limits<WT>::max)(),
                                            zero);
  }

  namespace detail {

    template <class VertexAndEdgeListGraph, class DistanceMatrix,
              class P, class T, class R, class Weight, 
              class VertexID>
    bool
    johnson_dispatch(VertexAndEdgeListGraph& g, 
                     DistanceMatrix& D,
                     const bgl_named_params<P, T, R>& params,
                     Weight w, VertexID id)
    {
      typedef typename property_traits<Weight>::value_type WT;
      
      return johnson_all_pairs_shortest_paths
        (g, D, id, w,
        choose_param(get_param(params, distance_compare_t()), 
          std::less<WT>()),
        choose_param(get_param(params, distance_combine_t()), 
          closed_plus<WT>()),
        choose_param(get_param(params, distance_inf_t()), 
          std::numeric_limits<WT>::max BOOST_PREVENT_MACRO_SUBSTITUTION()),
         choose_param(get_param(params, distance_zero_t()), WT()) );
    }

  } // namespace detail

  template <class VertexAndEdgeListGraph, class DistanceMatrix,
            class P, class T, class R>
  bool
  johnson_all_pairs_shortest_paths
    (VertexAndEdgeListGraph& g, 
     DistanceMatrix& D,
     const bgl_named_params<P, T, R>& params)
  {
    return detail::johnson_dispatch
      (g, D, params,
       choose_const_pmap(get_param(params, edge_weight), g, edge_weight),
       choose_const_pmap(get_param(params, vertex_index), g, vertex_index)
       );
  }

  template <class VertexAndEdgeListGraph, class DistanceMatrix>
  bool
  johnson_all_pairs_shortest_paths
    (VertexAndEdgeListGraph& g, DistanceMatrix& D)
  {
    bgl_named_params<int,int> params(1);
    return detail::johnson_dispatch
      (g, D, params, get(edge_weight, g), get(vertex_index, g));
  }

} // namespace boost

#endif // BOOST_GRAPH_JOHNSON_HPP


