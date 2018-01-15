/**
 * @file searches.cc
 *
 * Implement the BGL dfs and bfs wrappers.
 */

/*
 * David Gleich
 * 19 April 2006
 */

/*
 * 18 April 2007
 * Added src/dst pairs to all algorithms to stop searches early.
 * Corrected small documentation bugs.
 *
 * 9 July 2007
 * Switched to simple_csr_matrix graph type
 *
 * 29 August 2007
 * Switched to make_iterator_property_map in astar_search* functions
 * to fix compile bugs on g++-4.1
 */

#include "include/matlab_bgl.h"

#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <yasmic/iterator_utility.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/astar_search.hpp>
#include <utility>

#include "visitor_macros.hpp"
#include "stop_visitors.hpp"

struct stop_bfs {}; // stop bfs exception

/**
 * Wrap a boost graph library call to bfs.
 *
 * the ja and ia arrays specify the connectivity of the underlying graph,
 * ia is a length (nverts+1) array with the indices in ja that start the
 * nonzeros in each row.  ja is a length (ia(nverts)) array with the
 * columns of the connectivity.
 *
 * if d, dt, or pred is NULL, then that parameter is not computed.
 *   (currently not implemented)
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param src the source vertex for the search
 * @param dst a target vertex unless src=dst
 * @param d the distance array
 * @param dt the discover time array
 * @param pred the predecessor array
 * @return an error code if possible
 */
int breadth_first_search(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    int* d, int* dt, mbglIndex* pred /* output data: distance, discover time, predecessor */
    )
{
    using namespace yasmic;
    using namespace boost;
    using namespace std;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    if (!d || !dt || !pred)
    {
        return (-1);
    }

    // set the predecessor to itself
    for (mbglIndex i = 0; i < nverts; i++)
    {
        pred[i] = i;
    }

    int time = 0;

    if (dst == nverts)
    {
        // this case means we don't stop the search
        boost::breadth_first_search(g, src,
                    boost::visitor(make_bfs_visitor(
                        make_pair(record_distances(d, on_tree_edge()),
                        make_pair(stamp_times(dt, time, on_discover_vertex()),
                                  record_predecessors(pred, on_tree_edge()))))));
    }
    else
    {
        // this case means we do stop the search
        try {
            boost::breadth_first_search(g, src,
                        boost::visitor(make_bfs_visitor(
                            make_pair(record_distances(d, on_tree_edge()),
                            make_pair(stamp_times(dt, time, on_discover_vertex()),
                            make_pair(record_predecessors(pred, on_tree_edge()),
                                stop_search_on_vertex_target(dst, stop_bfs(), on_discover_vertex())))))));
        } catch (stop_bfs) {}
    }

    return (0);
}

template <class Graph>
struct c_bfs_visitor
{
    bfs_visitor_funcs_t *vis;

    VISITOR_VERTEX_FUNC(initialize_vertex, stop_bfs)
    VISITOR_VERTEX_FUNC(examine_vertex, stop_bfs)
    VISITOR_VERTEX_FUNC(discover_vertex, stop_bfs)
    VISITOR_VERTEX_FUNC(finish_vertex, stop_bfs)

    VISITOR_EDGE_FUNC(examine_edge, stop_bfs)
    VISITOR_EDGE_FUNC(tree_edge, stop_bfs)
    VISITOR_EDGE_FUNC(non_tree_edge, stop_bfs)
    VISITOR_EDGE_FUNC(gray_target, stop_bfs)
    VISITOR_EDGE_FUNC(black_target, stop_bfs)
};


int breadth_first_search_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, /* problem data */
    bfs_visitor_funcs_t vis /* visitor */
    )
{
    /*
     * History
     *
     * 2007-04-17
     * Added try-catch for stop_bfs exception
     */
    using namespace yasmic;
    using namespace boost;
    using namespace std;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    c_bfs_visitor<crs_graph> visitor;
    visitor.vis = &vis;

    try {
        boost::breadth_first_search(g, src,
            boost::visitor(visitor));
    }
    catch (stop_bfs)
    {
    }

    return (0);
}


struct stop_dfs {}; // stop dfs exception

/**
 * Wrap a boost graph library call to dfs.
 *
 * the ja and ia arrays specify the connectivity of the underlying graph,
 * ia is a length (nverts+1) array with the indices in ja that start the
 * nonzeros in each row.  ja is a length (ia(nverts)) array with the
 * columns of the connectivity.
 *
 * if dt, ft, or pred is NULL, then that parameter is not computed.
 *   (currently not implemented)
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param src the source vertex for the flow
 * @param full if full is non-zero, then compute the full dfs over all
 *    vertices, not just the connected component.
 * @param dst the target vertex if n < nverts and src != dst which
 *    stops the dfs if it is reached
 * @param d
 * @param dt the discover time array
 * @param ft the finish time array
 * @param pred the predecessor array
 * @return an error code if possible
 */
int depth_first_search(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, int full, mbglIndex dst, /* problem data */
    int* d, int* dt, int *ft, mbglIndex *pred /* output data: discover time, finish time, predecessor */
    )
{
    using namespace yasmic;
    using namespace boost;
    using namespace std;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    if (!d || !dt || !ft || !pred)
    {
        return (-1);
    }

    // set the predecessor to itself
    for (mbglIndex i = 0; i < nverts; i++)
    {
        pred[i] = i;
    }

    int time=0;

    if (!full)
    {
        std::vector<default_color_type> color_vec(num_vertices(g));
        default_color_type c = white_color; // avoid warning about un-init

        //
        // call visit because they don't want the full dfs
        //

        if (dst == nverts) {
            // they don't want to stop early
            boost::depth_first_visit(g, src,
                        make_dfs_visitor(
                            make_pair(record_distances(d, on_tree_edge()),
                            make_pair(stamp_times(dt, time, on_discover_vertex()),
                            make_pair(stamp_times(ft, time, on_finish_vertex()),
                                      record_predecessors(pred, on_tree_edge()))))),
                make_iterator_property_map(color_vec.begin(),
                    get(vertex_index,g)));
        } else {
            // they want to try and stop early
            try {
                boost::depth_first_visit(g, src,
                            make_dfs_visitor(
                                make_pair(record_distances(d, on_tree_edge()),
                                make_pair(stamp_times(dt, time, on_discover_vertex()),
                                make_pair(stamp_times(ft, time, on_finish_vertex()),
                                make_pair(record_predecessors(pred, on_tree_edge()),
                                    stop_search_on_vertex_target(dst, stop_dfs(), on_discover_vertex())))))),
                    make_iterator_property_map(color_vec.begin(),
                        get(vertex_index,g)));
            } catch (stop_dfs) {}
        }
    }
    else
    {
        // call search because they want the full dfs
        boost::depth_first_search(g,
                    boost::visitor(make_dfs_visitor(
                        make_pair(record_distances(d, on_tree_edge()),
                        make_pair(stamp_times(dt, time, on_discover_vertex()),
                                  record_predecessors(pred, on_tree_edge()))))));
    }
    return (0);
}

template <class Graph>
struct c_dfs_visitor
{
    dfs_visitor_funcs_t *vis;

    VISITOR_VERTEX_FUNC(initialize_vertex, stop_dfs)
    VISITOR_VERTEX_FUNC(start_vertex, stop_dfs)
    VISITOR_VERTEX_FUNC(discover_vertex, stop_dfs)
    VISITOR_VERTEX_FUNC(finish_vertex, stop_dfs)

    VISITOR_EDGE_FUNC(examine_edge, stop_dfs)
    VISITOR_EDGE_FUNC(tree_edge, stop_dfs)
    VISITOR_EDGE_FUNC(back_edge, stop_dfs)
    VISITOR_EDGE_FUNC(forward_or_cross_edge, stop_dfs)
};

int depth_first_search_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, int full, /* problem data */
    dfs_visitor_funcs_t vis)
{
    using namespace yasmic;
    using namespace boost;
    using namespace std;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    c_dfs_visitor<crs_graph> visitor;
    visitor.vis = &vis;

    try
    {
        if (!full)
        {
            std::vector<default_color_type> color_vec(num_vertices(g));
            default_color_type c = white_color; // avoid warning about un-init

            // call visit because they don't want the full dfs
            boost::depth_first_visit(g, src,
                visitor,
                make_iterator_property_map(color_vec.begin(),
                    get(vertex_index,g)));
        }
        else
        {
            // call search because they want the full dfs
            boost::depth_first_search(g,
                        boost::visitor(visitor)
                        );
        }
    }
    catch (stop_dfs)
    {
    }

    return (0);
}

struct stop_astar {}; // stop astar exception

template <class Graph>
struct c_astar_visitor
{
    astar_visitor_funcs_t *vis;

    VISITOR_VERTEX_FUNC(initialize_vertex, stop_astar)
    VISITOR_VERTEX_FUNC(examine_vertex, stop_astar)
    VISITOR_VERTEX_FUNC(discover_vertex, stop_astar)
    VISITOR_VERTEX_FUNC(finish_vertex, stop_astar)

    VISITOR_EDGE_FUNC(examine_edge, stop_astar)
    VISITOR_EDGE_FUNC(edge_relaxed, stop_astar)
    VISITOR_EDGE_FUNC(edge_not_relaxed, stop_astar)
    VISITOR_EDGE_FUNC(black_target, stop_astar)
};

template <class Graph, class CostType>
class astar_heuristic_data : public std::unary_function<
    typename boost::graph_traits<Graph>::vertex_descriptor, CostType>
{
private:
    CostType *_data;

public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    astar_heuristic_data(CostType *data) : _data(data) {}
    CostType operator()(Vertex u) { return _data[u]; }
};

template <class Graph, class CostType>
class astar_heuristic_func : public std::unary_function<
    typename boost::graph_traits<Graph>::vertex_descriptor, CostType>
{
private:
    void *_pdata;
    CostType (*_func)(void *pdata, mbglIndex u);

public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    astar_heuristic_func(CostType (*hfunc)(void* pdata, mbglIndex u), void* pdata) : _func(hfunc), _pdata(pdata) {}
    CostType operator()(Vertex u) { return _func(_pdata, u); }
};


int astar_search(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* start, end vertex */
    double *d, mbglIndex *pred, double *f, /* output */
    double *h /* heuristic function value for all vertices */, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    // create the graph g
    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    astar_heuristic_data<crs_weighted_graph,double> hdata(h);

    if (dst == nverts)
    {
        astar_search(g, src, hdata,
            distance_inf(dinf).
            predecessor_map(make_iterator_property_map(pred, get(vertex_index,g))).
            rank_map(make_iterator_property_map(f, get(vertex_index,g))).
            distance_map(make_iterator_property_map(d, get(vertex_index,g))));
    } else {
        try {
            astar_search(g, src, hdata,
                distance_inf(dinf).
                predecessor_map(make_iterator_property_map(pred, get(vertex_index,g))).
                rank_map(make_iterator_property_map(f, get(vertex_index,g))).
                distance_map(make_iterator_property_map(d, get(vertex_index,g))).
                visitor(make_astar_visitor(
                    stop_search_on_vertex_target(dst, stop_astar(), on_discover_vertex()))));
        } catch (stop_astar) {}
    }

    return (0);
}

int astar_search_hfunc(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* start vertex */
    double *d, mbglIndex *pred, double *f,
    double (*hfunc)(void* pdata, mbglIndex u) /* heuristic function */, void* pdata, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    astar_heuristic_func<crs_weighted_graph,double> h(hfunc, pdata);

    if (dst == nverts)
    {
        astar_search(g, src, h,
            distance_inf(dinf).
            predecessor_map(make_iterator_property_map(pred, get(vertex_index,g))).
            rank_map(make_iterator_property_map(f, get(vertex_index,g))).
            distance_map(make_iterator_property_map(d, get(vertex_index,g))));
    } else {
        try {
            astar_search(g, src, h,
                distance_inf(dinf).
                predecessor_map(make_iterator_property_map(pred, get(vertex_index,g))).
                rank_map(make_iterator_property_map(f, get(vertex_index,g))).
                distance_map(make_iterator_property_map(d, get(vertex_index,g))).
                visitor(make_astar_visitor(
                    stop_search_on_vertex_target(dst, stop_astar(), on_discover_vertex()))));
        } catch (stop_astar) {}
    }

    return (0);
}

int astar_search_hfunc_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, /* start vertex */
    double *d, mbglIndex *pred, double *f,
    double (*hfunc)(void* pdata, mbglIndex u) /* heuristic function */, void* pdata, double dinf,
    astar_visitor_funcs_t vis)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    astar_heuristic_func<crs_weighted_graph,double> h(hfunc, pdata);

    c_astar_visitor<crs_weighted_graph> visitor_impl;
    visitor_impl.vis = &vis;

    try
    {
        astar_search(g, src, h,
            distance_inf(dinf).
            predecessor_map(make_iterator_property_map(pred, get(vertex_index,g))).
            rank_map(make_iterator_property_map(f, get(vertex_index,g))).
            distance_map(make_iterator_property_map(d, get(vertex_index,g))).
            visitor(visitor_impl));
    }
    catch (stop_astar)
    {
    }

    return (0);
}
