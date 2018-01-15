/** @file planar.cc
 * @author David F. Gleich
 * @date 2008-09-29
 * @copyright Stanford University, 2008
 * Planar graph algorithm wrappers
 */

/** History
 *  2008-09-29: Initial coding
 */

#include "include/matlab_bgl.h"

#include <yasmic/undir_simple_csr_matrix_as_graph.hpp>
#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <yasmic/iterator_utility.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
//#include <boost/graph/chrobak_payne_drawing.hpp>
#include <yasmic/boost_mod/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/is_kuratowski_subgraph.hpp>
#include <boost/graph/make_connected.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/simple_point.hpp>
//#include <boost/graph/is_straight_line_drawing.hpp>
#include <yasmic/boost_mod/is_straight_line_drawing.hpp>
#include <boost/graph/filtered_graph.hpp>

#include <vector>
#include <iostream>
#include <algorithm>

#include <math.h>

#include "libmbgl_util.hpp"

/** Copy an edge iterator to a pair of arrays
 */
template <typename Graph, typename Iterator>
mbglIndex copy_to_ij(Graph& g, Iterator oi,
                   Iterator oi_end, mbglIndex* i, mbglIndex* j)
{
  using namespace boost;

  mbglIndex ei;
  for (ei= 0; oi != oi_end; ++oi, ++ei) {
    typename graph_traits<Graph>::edge_descriptor e = *oi;
    i[ei] = source(e,g);
    j[ei] = target(e,g);
  }

  return (ei);
}

/** Translate the boost embedding information to the libmbgl embedding output
 */
template <typename Graph, typename PlanarEmbedding>
void copy_embedding(Graph& g, PlanarEmbedding e, mbglIndex *eip, mbglIndex *eie)
{
  using namespace boost;
  typename graph_traits<Graph>::vertex_iterator vi, viend;
  mbglIndex nedges= num_edges(g);
  mbglIndex curi= 0;
  mbglIndex nullv= graph_traits<Graph>::null_vertex();
  mbglIndex oldv= nullv;
  for (boost::tie(vi,viend)=vertices(g); vi!=viend; ++vi) {
    if (oldv!=nullv) { assert(*vi==oldv+1); }
    assert(*vi<num_vertices(g));
    eip[*vi] = curi;
    typename property_traits<PlanarEmbedding>::value_type::const_iterator
      ei=e[*vi].begin(), eiend=e[*vi].end();
    for (; ei!=eiend; ++ei) {
      assert(source(*ei,g) == *vi || target(*ei,g) == *vi);
      assert(curi<nedges);
      if (source(*ei,g) == *vi) {
        eie[curi] = target(*ei,g);
      } else if (target(*ei,g) == *vi) {
        eie[curi] = source(*ei,g);
      } else {
        assert(source(*ei,g) == *vi || target(*ei,g) == *vi);
      }
      ++curi;
    }
    oldv = *vi;
  }
  eip[num_vertices(g)] = curi;
}

/** Test if a graph is planar, compute a planar embedding, or get a Kuratowski
 * subgraph
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param is_planar set to 0 if the graph is not planar, else the graph is
 *   planar
 * @param i the source vertex of any edge in the Kuratowski subgraph,
 *   length max(3*nverts-6,6)
 * @param j the dest vertex of any edge in the Kuratowski subgraph,
 *   length max(3*nverts-6,6)
 * @param nedges set to the number of edges in i and j actually used
 * @param eip the embedding information edge pointer to eie, length nverts+1
 *   eip[eie[v]] to eip[eie[v+1]] gives the embedding order for vertex v
 * @param eie the embedding information destination list, length ia[nverts]
 *   see eip for a description.
 */
int boyer_myrvold_planarity_test(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int *is_planar,
    mbglIndex *i, mbglIndex *j, mbglIndex* nedges, /* kuratowski subgraph output */
    mbglIndex *eip, mbglIndex *eie)
{
  using namespace yasmic;
  using namespace boost;

  typedef simple_csr_matrix<mbglIndex,double> crs_graph;
  crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

  //Define the storage type for the planar embedding
  typedef std::vector< std::vector< graph_traits<crs_graph>::edge_descriptor > >
    embedding_storage_t;
  typedef boost::iterator_property_map
    < embedding_storage_t::iterator,
      property_map<crs_graph, vertex_index_t>::type
    >
    embedding_t;


  bool planar=false;
  if ((i == NULL || j == NULL || nedges == NULL) &&
      (eip == NULL || eip == NULL)) {
    // just test for a planar graph
    planar= boyer_myrvold_planarity_test(g);
  } else if (eip == NULL || eie == NULL ) {
    // just get the kuratowski subgraph
    std::vector<graph_traits<crs_graph>::edge_descriptor>
      kuratowski_edges;
    planar= boyer_myrvold_planarity_test(
        boyer_myrvold_params::graph = g,
        boyer_myrvold_params::kuratowski_subgraph =
          std::back_inserter(kuratowski_edges)
        );
    if (planar == false) {
      *nedges =
        copy_to_ij(g, kuratowski_edges.begin(), kuratowski_edges.end(), i, j);
    } else { *nedges = 0; }
  } else if (i == NULL || j == NULL || nedges == NULL) {
    embedding_storage_t embedding_storage(num_vertices(g));
    embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));
    planar= boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                   boyer_myrvold_params::embedding = embedding
                                   );
    if (planar == true) {
      copy_embedding(g, embedding, eip, eie);
    } else {
      memset(eip, 0, sizeof(mbglIndex)*(nverts+1));
    }
  } else {
    // get the kuratowski subgraph and the planar embedding
    std::vector<graph_traits<crs_graph>::edge_descriptor>
        kuratowski_edges;
    embedding_storage_t embedding_storage(num_vertices(g));
    embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));
    planar= boyer_myrvold_planarity_test(
        boyer_myrvold_params::graph = g,
        boyer_myrvold_params::embedding = embedding,
        boyer_myrvold_params::kuratowski_subgraph =
                 std::back_inserter(kuratowski_edges)
        );
    if (planar == false) {
      *nedges =
        copy_to_ij(g, kuratowski_edges.begin(), kuratowski_edges.end(), i, j);
    } else { *nedges = 0; }
    if (planar == true) {
      copy_embedding(g, embedding, eip, eie);
    } else {
      memset(eip, 0, sizeof(mbglIndex)*(nverts+1));
    }
  }
  if (is_planar) {
    if (planar) { *is_planar = 1; }
    else { *is_planar = 0; }
  }
  return (0);
}

/** An edge filter that only returns the edges in an upper-triangular part
 */
template <typename Graph>
struct upper_triangle_edge_filter {
  upper_triangle_edge_filter() : _g(NULL) {}
  upper_triangle_edge_filter(const Graph& g_) : _g(&g_) {}
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return (boost::source(e,*_g) < boost::target(e,*_g));
  }
  const Graph* _g;
};


/** Test is a graph is Kuratowski (i.e. contracts to K_3,3 or K_5)
 * This function only uses the upper triangular set of edges.  It will
 * silently report incorrect ouput if there are too many parallel edges,
 * so beware.
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param is_ksubgraph set to 1 if the graph is Kuratowski
 */
int is_kuratowski_subgraph(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int *is_ksubgraph)
{
  using namespace boost;
  using namespace yasmic;

  typedef simple_csr_matrix<mbglIndex,double> crs_graph;
  typedef upper_triangle_edge_filter<crs_graph> filter;
  typedef filtered_graph<crs_graph, filter> fgraph;
  assert(is_ksubgraph);
  crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);
  filter f(g);
  fgraph fg(g,f);
  graph_traits<fgraph>::edge_iterator ei, eiend;
  boost::tie(ei,eiend)= edges(fg);
  bool is_k= is_kuratowski_subgraph(g, ei, eiend, get(vertex_index,g));
  if (is_k) { *is_ksubgraph = 1; } else { *is_ksubgraph = 0; }
  return (0);
}

/** Test if a given set of positions is a straight line drawing of the graph.
 * This function uses a bucket-sort, and so large positions may cause out
 * of memory errors.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param X The set of positions.  These are rounded to size_t variables and so
 *   all the elements of X must be positive.
 * @param is_sldrawing set to 1 if the positions are a straight line drawing,
 *   otherwise set to 0
 */
int is_straight_line_drawing(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    double *X, int *is_sldrawing)
{
  using namespace boost;
  using namespace yasmic;

  typedef simple_csr_matrix<mbglIndex,double> crs_graph;
  crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);
  assert(is_sldrawing);
  // copy the layout from positions
  std::vector<simple_point<size_t> > position_vec(nverts);
  mbglIndex n = num_vertices(g);
  for (mbglIndex i = 0; i<n; i++) {
    if (X[i+0*n] < 0 || X[i+1*n] < 0) {
      return -11;
    }
    position_vec[i].x = (size_t)floor(X[i+0*n]);
    position_vec[i].y = (size_t)floor(X[i+1*n]);
  }
  bool is_sl= is_straight_line_drawing(g,
      make_iterator_property_map(position_vec.begin(),get(vertex_index,g)),
      get(vertex_index,g));
  if (is_sl) { *is_sldrawing = 1; } else { *is_sldrawing = 0; }
  return (0);
}

/** Copy a crs_graph to a mutable boost graph
 */
template <typename Graph>
void copy_crs_to_graph(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    Graph& g)
{
  using namespace boost;
  using namespace yasmic;
  typedef simple_csr_matrix<mbglIndex,double> crs_graph;
  crs_graph g_in(nverts, nverts, ia[nverts], ia, ja, NULL);
  //typename property_map<Graph, vertex_index_t>::type vi(get(vertex_index, g));
  //for (mbglIndex i=0; i<nverts; i++) { vi[i]=i; }
  typename graph_traits<crs_graph>::vertex_iterator vi, vi_end;
  typename graph_traits<crs_graph>::out_edge_iterator ei, ei_end;
  for (tie(vi, vi_end) = vertices(g_in); vi != vi_end; ++vi) {
    for (tie(ei, ei_end) = out_edges(*vi, g_in); ei != ei_end; ++ei) {
      if (source(*ei,g_in) < target(*ei,g_in)) {
        add_edge(source(*ei,g_in),target(*ei,g_in),g);
      }
    }
  }
}

/** Add an edge and record all edges we added
 * EdgeSrcOutIterator - the type of the out iterator for the source vertex
 *   in the new edge pair
 * EdgeDstOutIterator - the type of the out iterator for the destination vertex
 *   in the new edge pair
 * TODO Make this function update the edge_index parameter too!  This will avoid
 * O(E) work in a few cases.
 */
template <typename EdgeSrcOutIterator, typename EdgeDstOutIterator>
struct record_add_edge_visitor {
  EdgeSrcOutIterator soi, soi_end;
  EdgeDstOutIterator doi, doi_end;
  mbglIndex* count; // incremented for each additional edge

  record_add_edge_visitor(EdgeSrcOutIterator soi_, EdgeSrcOutIterator soi_end_,
      EdgeDstOutIterator doi_, EdgeDstOutIterator doi_end_, mbglIndex* count_)
  : soi(soi_), soi_end(soi_end_), doi(doi_), doi_end(doi_end_), count(count_)
  {}

  template <typename Graph, typename Vertex>
  void visit_vertex_pair(Vertex u, Vertex v, Graph& g) {
    add_edge(u,v,g);
    assert(soi != soi_end);
    assert(doi != doi_end);
    *soi = u;
    *doi = v;
    ++doi;
    ++soi;
    (*count)++;
  }
};

/** Helper to abstract the details of the triangulation */
template <typename Graph, typename RecordAddEdgeVisitor>
int triangulate_bgl_graph(
  Graph& g, int make_conn, int make_biconnected, int make_maximal,
  RecordAddEdgeVisitor add_edge_visitor)
{
  using namespace boost;
  typedef std::vector< std::vector<
            typename graph_traits<Graph>::edge_descriptor > >
  embedding_storage_t;
  typedef boost::iterator_property_map
    < typename embedding_storage_t::iterator,
      typename property_map<Graph, vertex_index_t>::type
    >
    embedding_t;

  if (make_conn) {
    make_connected(g, get(vertex_index,g), add_edge_visitor);
  }
  if (make_biconnected || make_maximal) {
    typename property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    typename graph_traits<Graph>::edges_size_type edge_count = 0;
    typename graph_traits<Graph>::edge_iterator ei, ei_end;

    // compute a planar embedding
    embedding_storage_t embedding_storage(num_vertices(g));
    embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));

    if (make_biconnected) {
      // compute the edge index
      edge_count = 0;
      for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        put(e_index, *ei, edge_count++);
      }
      // compute a planar embedding
      if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                       boyer_myrvold_params::embedding =
                                           embedding
                                       )
          ) {
        make_biconnected_planar(g, embedding, get(edge_index,g),
            add_edge_visitor);
      } else {
        return 1;
      }
    }
    if (make_maximal) {
      // compute the edge index
      edge_count = 0;
      for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        put(e_index, *ei, edge_count++);
      }
      // compute a planar embedding
      if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                       boyer_myrvold_params::embedding =
                                           embedding
                                       )
          ) {
        make_maximal_planar(g, embedding, get(vertex_index,g),
            get(edge_index,g), add_edge_visitor);
      } else {
        return 1;
      }
    }
  }
  return 0;
}


/** Compute extra edges that are needed for a straight line embedding.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param make_connected if set (!=0), compute edges to make a single
 *   connected component
 * @param make_biconnected If set (!=0), compute edges to make a single
 *   biconnected component (requires a planar input graph).  If only this
 *   parameter is set, then we assume the graph is already connected.
 * @param make_maximal If set (!=0), compute edges to make the maximal
 *   planar graph (requires a planar input graph). If only this
 *   parameter is set, then we assume the graph is already biconnected.
 * @param i source of all edges added in the triangulation,
 *   length max(3*nverts-6,6)
 * @param j dest of all edges added in the triangulation,
 *   length max(3*nverts-6,6)
 * @param nedges the number of edges in i and j actually used
 */
int triangulate_graph(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int make_connected, int make_biconnected, int make_maximal,
    mbglIndex *i, mbglIndex *j, mbglIndex* nedges /* extra edges */)
{
  using namespace boost;
  typedef adjacency_list
        < vecS, vecS, undirectedS, no_property,
        property<edge_index_t, mbglIndex> > graph;
  assert(i != NULL); assert(j != NULL); assert(nedges != NULL);
  *nedges = 0;
  size_t endedges= nverts>4 ? 3*nverts-6 : 6;
  record_add_edge_visitor<mbglIndex*,mbglIndex*> add_edge_visitor(
      i, i+endedges, j, j+endedges, nedges);
  graph g(nverts);
  copy_crs_to_graph(nverts, ja, ia, g);
  return triangulate_bgl_graph(
      g, make_connected, make_biconnected, make_maximal, add_edge_visitor);
}

/** Abstract some details between the csr_graph and the adjacency_list graph
 * for the straight line drawing
 */
template <typename Graph>
int chrobak_payne_straight_line_drawing_on_maximal_graph(const Graph& g,
    int just_ordering,
    mbglIndex *p, /* ordering permutation */
    mbglDegreeType *X)
{
  using namespace boost;
  typedef std::vector< std::vector<
            typename graph_traits<Graph>::edge_descriptor > >
    embedding_storage_t;
  typedef boost::iterator_property_map
    < typename embedding_storage_t::iterator,
      typename property_map<Graph, vertex_index_t>::type
    >
    embedding_t;

  embedding_storage_t embedding_storage(num_vertices(g));
  embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));

  if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                   boyer_myrvold_params::embedding =
                                      embedding
                                  )
     ) {
   // do nothing, the graph is planar
  } else {
   return 1;
  }
  // TODO convert this to use p directly
  mbglIndex n = num_vertices(g);
  assert(n == 0 || p);
  std::vector< typename graph_traits<Graph>::vertex_descriptor > ordering;
  planar_canonical_ordering(g, embedding,
     std::back_inserter(ordering), get(vertex_index,g));

  if (just_ordering == 0) {
    typedef simple_point<size_t> coord_t;
    typedef std::vector< coord_t > straight_line_drawing_storage_t;
    typedef boost::iterator_property_map
       < straight_line_drawing_storage_t::iterator,
         typename property_map<Graph, vertex_index_t>::type
       >
       straight_line_drawing_t;

    assert(n == 0 || X);
    // compute the straight line embedding
    straight_line_drawing_storage_t straight_line_drawing_storage
     (num_vertices(g));
    straight_line_drawing_t straight_line_drawing
     (straight_line_drawing_storage.begin(),
      get(vertex_index,g)
      );
    if (n == 2) {
      straight_line_drawing_storage[0].x = 0;
      straight_line_drawing_storage[0].y = 0;
      straight_line_drawing_storage[1].x = 1;
      straight_line_drawing_storage[1].y = 0;
    } else {
      chrobak_payne_straight_line_drawing(g,
                                         embedding,
                                         ordering.begin(),
                                         ordering.end(),
                                         straight_line_drawing
                                         );
    }
    // copy all the data
    for (mbglIndex i= 0; i<n; i++) {
      X[i+0*n] = straight_line_drawing_storage[i].x;
      X[i+1*n] = straight_line_drawing_storage[i].y;
    }
  }
  assert((ordering.end() - ordering.begin()) == (ptrdiff_t)n);
  for (mbglIndex i= 0; i<n; i++) {
    p[i] = ordering[i];
  }
  return 0;
}

/** Compute a straight line embedding of a graph or a canonical planar
 * ordering
 *
 * This function must copy the graph unless is_maximal is set.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param just_ordering if set (!=0), then just compute the ordering permutation
 * @param is_maximal if set (!=0), assume the input is a maximal planar graph
 * @param i source edges for the maximal planar graph, length 3n-6
 * @param j dest edges for the maximal planar graph, length 3n-6
 * @param nedges the number of entries in i,j actually used, length 1
 * @param p the ordering permutation, length nverts
 * @param X the positions, length 2*nverts
 */
int chrobak_payne_straight_line_drawing(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int just_ordering, int is_maximal,
    mbglIndex *i, mbglIndex *j, mbglIndex* nedges, /* extra edges */
    mbglIndex *p, /* ordering permutation */
    mbglDegreeType *X)
{
  using namespace boost;
  using namespace yasmic;
  if (is_maximal) {
    // trust the user that the graph is maximal, this means we
    // can skip the copy
    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);
    if (nedges) { *nedges = 0; } // already maximal
    return chrobak_payne_straight_line_drawing_on_maximal_graph(g,
        just_ordering, p, X);
  } else {
    typedef adjacency_list
      < vecS, vecS, undirectedS, no_property,
      property<edge_index_t, mbglIndex> > graph;
    graph g(nverts);
    copy_crs_to_graph(nverts, ja, ia, g);
    int rval=0;
    if (i != NULL && j != NULL && nedges != NULL) {
      *nedges = 0;
      size_t endedges= nverts>4 ? 3*nverts-6 : 6;
      record_add_edge_visitor<mbglIndex*,mbglIndex*> add_edge_visitor(
            i, i+endedges, j, j+endedges, nedges);
      rval= triangulate_bgl_graph(g, 1, 1, 1, add_edge_visitor);
    } else {
      mbglIndex cval1= 0;
      mbglIndex extra_edges= 0;
      typedef yasmic::constant_iterator<mbglIndex> dummy_iterator;
      dummy_iterator di(&cval1);
      record_add_edge_visitor<dummy_iterator,dummy_iterator> add_edge_visitor(
            di, di, di, di, &extra_edges);
      rval= triangulate_bgl_graph(g, 1, 1, 1, add_edge_visitor);
    }
    if (rval != 0) { return rval; } // return on error
    return chrobak_payne_straight_line_drawing_on_maximal_graph(g,
            just_ordering, p, X);
  }
}





