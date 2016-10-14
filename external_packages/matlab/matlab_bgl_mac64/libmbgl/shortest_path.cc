/**
 * @file shortest_path.cc
 *
 * Implement the BGL shortest path algorithm wrappers.
 */

/*
 * David Gleich
 * 19 April 2006
 */

/*
 * 18 April 2007
 * Added src/dst vertex pairs to all the calls to allow partial searches.
 * Corrected small documentation bugs.
 *
 * 9 July 2007
 * Switched to simple_csr_matrix graph type
 */

#include "include/matlab_bgl.h"

#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <yasmic/iterator_utility.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <yasmic/boost_mod/bellman_ford_shortest_paths.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <yasmic/boost_mod/floyd_warshall_shortest.hpp>

#include "visitor_macros.hpp"
#include "stop_visitors.hpp"
#include "libmbgl_util.hpp"

struct stop_dijkstra {}; // stop dijkstra exception

int dijkstra_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    if (dst == nverts) {
        dijkstra_shortest_paths(g, src, distance_inf(dinf).predecessor_map(pred).distance_map(d));
    } else {
        try {
            dijkstra_shortest_paths(g, src,
                distance_inf(dinf).predecessor_map(pred).distance_map(d).
                visitor(make_dijkstra_visitor(
                            stop_search_on_vertex_target(dst, stop_dijkstra(), on_discover_vertex()))));
        } catch (stop_dijkstra) {}
    }

    return (0);
}

template <class Graph>
struct c_dijkstra_visitor
{
    dijkstra_visitor_funcs_t *vis;

    VISITOR_VERTEX_FUNC(initialize_vertex, stop_dijkstra)
    VISITOR_VERTEX_FUNC(examine_vertex, stop_dijkstra)
    VISITOR_VERTEX_FUNC(discover_vertex, stop_dijkstra)
    VISITOR_VERTEX_FUNC(finish_vertex, stop_dijkstra)

    VISITOR_EDGE_FUNC(examine_edge, stop_dijkstra)
    VISITOR_EDGE_FUNC(edge_relaxed, stop_dijkstra)
    VISITOR_EDGE_FUNC(edge_not_relaxed, stop_dijkstra)

};

int dijkstra_sp_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, /* problem data */
    double* d, mbglIndex *pred,
    double dinf, dijkstra_visitor_funcs_t vis)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    c_dijkstra_visitor<crs_weighted_graph> visitor_impl;
    visitor_impl.vis = &vis;

    try
    {
        dijkstra_shortest_paths(g, src, distance_inf(dinf).predecessor_map(pred).distance_map(d).visitor(visitor_impl));
    }
    catch (stop_dijkstra)
    {
    }

    return (0);
}

struct stop_bellman_ford {}; // stop bellman_ford exception

int bellman_ford_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    if (dst == nverts) {
        bellman_ford_shortest_paths(g,
            root_vertex(src).distance_inf(dinf).predecessor_map(pred).distance_map(d));
    } else {
        try {
            bellman_ford_shortest_paths(g,
                root_vertex(src).distance_inf(dinf).predecessor_map(pred).distance_map(d).
                visitor(make_bellman_visitor(
                            stop_search_on_vertex_target(dst, stop_bellman_ford(), on_discover_vertex()))));
        } catch (stop_bellman_ford) {}
    }

    return (0);
}

template <class Graph>
struct c_bellman_ford_visitor
{
    bellman_ford_visitor_funcs_t *vis;

    VISITOR_VERTEX_FUNC(initialize_vertex, stop_bellman_ford)

    VISITOR_EDGE_FUNC(examine_edge, stop_bellman_ford)
    VISITOR_EDGE_FUNC(edge_relaxed, stop_bellman_ford)
    VISITOR_EDGE_FUNC(edge_not_relaxed, stop_bellman_ford)
    VISITOR_EDGE_FUNC(edge_minimized, stop_bellman_ford)
    VISITOR_EDGE_FUNC(edge_not_minimized, stop_bellman_ford)
};

int bellman_ford_sp_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src,  /* problem data */
    double* d, mbglIndex *pred, double dinf,
    bellman_ford_visitor_funcs_t vis)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    c_bellman_ford_visitor<crs_weighted_graph> visitor_impl;
    visitor_impl.vis = &vis;

    try
    {
        bellman_ford_shortest_paths(g,
            root_vertex(src).distance_inf(dinf).predecessor_map(pred).distance_map(d).visitor(visitor_impl));
    }
    catch (stop_bellman_ford)
    {
    }

    return (0);
}

struct stop_dag {}; // stop dag exception

int dag_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    if (dst == nverts) {
        dag_shortest_paths(g, src,
            distance_inf(dinf).predecessor_map(pred).distance_map(d));
    } else {
        try {
            dag_shortest_paths(g, src,
                distance_inf(dinf).predecessor_map(pred).distance_map(d).
                visitor(make_dijkstra_visitor(
                    stop_search_on_vertex_target(dst, stop_dag(), on_discover_vertex()))));
        } catch (stop_dag) {}
    }

    return (0);
}

int johnson_all_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double* D, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    row_matrix<double> Dmat(D,nverts,nverts);

    bool rval = johnson_all_pairs_shortest_paths(g, Dmat,
		distance_inf(dinf).distance_combine(std::plus<double>()));

	if (rval == true)
	{
		return (0);
	}

	// else, there was an error
	return (-1);
}

int floyd_warshall_all_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double* D, double dinf,
    mbglIndex* pred)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    row_matrix<double> Dmat(D,nverts,nverts);
    bool rval = false;
    if (!pred) {
        rval = floyd_warshall_all_pairs_shortest_paths(g,
            Dmat, distance_inf(dinf).distance_combine(std::plus<double>()));
    } else {
        row_matrix<mbglIndex> Pmat(pred,nverts,nverts);
        //rval = floyd_warshall_all_pairs_shortest_paths(g,
        //    Dmat, distance_inf(dinf).distance_combine(std::plus<double>()).predecessor_map(Pmat));
        // making this call is ridiculous, but otherwise, it won't pick up the right type!
        rval = floyd_warshall_all_pairs_shortest_paths(g, Dmat, Pmat, get(edge_weight,g),
            std::less<double>(),std::plus<double>(),dinf,0.0);
    }

	if (rval == true)
	{
		return (0);
	}

	// else, there was an error
	return (-1);
}

