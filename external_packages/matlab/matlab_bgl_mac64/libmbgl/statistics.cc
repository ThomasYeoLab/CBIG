/** @file statistics.cc
 * @author David F. Gleich
 * @date 2006-04-21
 * @copyright Stanford University, 2006-2008
 * Graph statistics wrappers
 */

/** History
 *  2007-04-18: Implemented topological order function
 *  2007-07-05: Implemented matchings
 *  2007-07-05: Implemented core_numbers
 *  2007-07-11: Implemented directed and weighted clustering coefficients
 *  2007-07-12: Implemented dominator tree
 */

#include "include/matlab_bgl.h"

#include <yasmic/simple_csr_matrix_as_graph.hpp>
#include <yasmic/simple_row_and_column_matrix_as_graph.hpp>
#include <yasmic/iterator_utility.hpp>

#include <vector>

#include <boost/graph/iteration_macros.hpp>
//#include <boost/graph/betweenness_centrality.hpp>
#include <yasmic/boost_mod/betweenness_centrality.hpp>
#include <boost/graph/topological_sort.hpp>

#include <boost/iterator/reverse_iterator.hpp>

#include <boost/graph/max_cardinality_matching.hpp>
#include <yasmic/boost_mod/core_numbers.hpp>
#include <boost/graph/dominator_tree.hpp>

#include <iostream>

#include <math.h>

template <class Vertex, class IndMap>
struct in_indicator_pred
	: public std::unary_function<Vertex, bool>
{
	in_indicator_pred(IndMap indmap, Vertex src)
		: i(indmap), u(src)
	{}

	IndMap i;
    Vertex u;

	bool operator() (const Vertex &v) const
	{
		return (i[v] > 0) && (u != v);
	}
};

/* Clustering Coefficients Code */
template <class Graph, class CCMap, class IndMap>
void cluster_coefficients(const Graph& g, CCMap cc, IndMap ind)
{
	using namespace boost;

	BGL_FORALL_VERTICES_T(v,g,Graph)
	{
		cc[v] = 0;
		ind[v] = 0;
	}


	BGL_FORALL_VERTICES_T(v,g,Graph)
	{
		BGL_FORALL_ADJ_T(v,w,g,Graph)
		{
			ind[w] = 1;
		}

        ind[v] = 0;

		typename property_traits<CCMap>::value_type cur_cc = 0;
		typename graph_traits<Graph>::degree_size_type d = out_degree(v, g);



		BGL_FORALL_ADJ_T(v,w,g,Graph)
		{
            // if we are adjacent to ourselves, skip the iteration
            if (v == w) { --d; continue; }

            in_indicator_pred<
    			typename graph_traits<Graph>::vertex_descriptor,
			    IndMap> p(ind,w);

			typename graph_traits<Graph>::adjacency_iterator ai, aiend;
			boost::tie(ai, aiend) = adjacent_vertices(w,g);

			// count if this is in the indicator predicate
			cur_cc += (int)count_if(ai, aiend, p);
		}

        if (d > 1)
        {
		    cc[v] = (double)cur_cc/(double)((d*(d-1)));
        }
        else
        {
            cc[v] = 0.0;
        }

		// reset the indicator
		BGL_FORALL_ADJ_T(v,w,g,Graph)
		{
			ind[w] = 0;
		}
	}
}

template <class Graph, class CCMap, class EdgeWeightMap, class IndMap>
void undirected_clustering_coefficients(const Graph& g, CCMap cc, EdgeWeightMap wm,
    IndMap ind)
{
    using namespace boost;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex;
    typename graph_traits<Graph>::vertex_iterator vi,vi_end;
    for (tie(vi,vi_end)=vertices(g);vi!=vi_end;++vi) {
        put(cc,*vi,0);
        put(ind,*vi,0);
    }
    // a lazy cache for the temp weights for each vertex, could be more
    // efficient using the maximum degree
    std::vector<typename property_traits<EdgeWeightMap>::value_type> cache(num_vertices(g));
    std::vector<typename graph_traits<Graph>::degree_size_type> degs;
    for (tie(vi,vi_end)=vertices(g);vi!=vi_end;++vi) {
        vertex v = *vi;
        typename graph_traits<Graph>::out_edge_iterator oi,oi_end;
        typename graph_traits<Graph>::out_edge_iterator oi2,oi2_end;
        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            // check to ignore self edges
            if (target(*oi,g) != v) {
                put(ind,target(*oi,g),1);
                cache[target(*oi,g)] = pow(get(wm,*oi),1.0/3.0);
            }
        }
        typename property_traits<CCMap>::value_type cur_cc = 0;
        typename graph_traits<Graph>::degree_size_type d = out_degree(v,g);
        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            vertex w = target(*oi,g);
            if (v==w) { --d; }
            for (tie(oi2,oi2_end)=out_edges(w,g);oi2!=oi2_end;++oi2) {
                if (target(*oi2,g) == w) { continue; }
                if (get(ind,target(*oi2,g))) {
                    // cache is cached with its power
                    cur_cc += pow(get(wm,*oi2),1.0/3.0)*pow(get(wm,*oi),1.0/3.0)*cache[target(*oi2,g)];
                }
            }
        }
        if (d > 1) {
            put(cc,v,cur_cc/((typename property_traits<CCMap>::value_type)(d*(d-1))));
        } else {
            put(cc,v,(typename property_traits<CCMap>::value_type)0);
        }
        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            put(ind,target(*oi,g),0);
            cache[target(*oi,g)] = 0;
        }
    }
}

// graph must be a bidirectional graph
template <class Graph, class CCMap, class EdgeWeightMap, class IndMap>
void directed_clustering_coefficients(const Graph& g,
    CCMap cc, EdgeWeightMap wm, IndMap ind)
{
    using namespace boost;
    typedef typename property_traits<CCMap>::value_type value_type;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex;
    typename graph_traits<Graph>::vertex_iterator vi,vi_end;
    for (tie(vi,vi_end)=vertices(g);vi!=vi_end;++vi) {
        put(cc,*vi,0);
        put(ind,*vi,0);
    }
    // a lazy cache for the temp weights for each vertex, could be more
    // efficient using the maximum degree
    std::vector<typename property_traits<EdgeWeightMap>::value_type> cache(num_vertices(g));
    std::vector<typename graph_traits<Graph>::degree_size_type> degs(num_vertices(g));
    typename graph_traits<Graph>::edge_iterator ei,ei_end;
    for (tie(ei,ei_end)=edges(g);ei!=ei_end;++ei) {
        vertex u = source(*ei,g);
        vertex v = target(*ei,g);
        if (u==v) { continue; }
        degs[u]++;
        degs[v]++;
    }
    for (tie(vi,vi_end)=vertices(g);vi!=vi_end;++vi) {
        vertex v = *vi;
        typename graph_traits<Graph>::out_edge_iterator oi,oi_end;
        typename graph_traits<Graph>::out_edge_iterator oi2,oi2_end;
        typename graph_traits<Graph>::in_edge_iterator ii, ii_end;
        for (tie(ii,ii_end)=in_edges(v,g);ii!=ii_end;++ii) {
            if (source(*ii,g) != v) {
                put(ind,source(*ii,g),1);
                cache[source(*ii,g)] += pow(get(wm,*ii),1.0/3.0);
            }
        }
        // we've precomputed the set of in-edges, so we know how
        // to get back to v to complete a triangle that ends
        // with an edge to v.
        typename graph_traits<Graph>::degree_size_type bilateral_edges = 0;
        value_type cur_cc_cyc=0, cur_cc_mid=0, cur_cc_in=0, cur_cc_out=0;
        // cycles are out-out-out.
        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            vertex w = target(*oi,g);
            if (v == w) { continue; }
            value_type cache_w = pow(get(wm,*oi),1.0/3.0);
            for (tie(oi2,oi2_end)=out_edges(w,g);oi2!=oi2_end;++oi2) {
                vertex u = target(*oi2,g);
                if (u == v) { ++bilateral_edges; continue; }
                if (u == w) { continue; }
                if (get(ind,u)) {
                    // cache is cached with its power
                    cur_cc_cyc += pow(get(wm,*oi2),1.0/3.0)*cache_w*cache[u];
                }
            }
        }
        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            vertex w = target(*oi,g);
            if (v == w) { continue; }
            value_type cache_w = pow(get(wm,*oi),1.0/3.0);
            for (tie(ii,ii_end)=in_edges(w,g);ii!=ii_end;++ii) {
                vertex u = source(*ii,g);
                if (u == w) { continue; }
                if (get(ind,u)) {
                    // cache is cached with its power
                    cur_cc_mid += pow(get(wm,*ii),1.0/3.0)*cache_w*cache[u];
                }
            }
        }
        for (tie(ii,ii_end)=in_edges(v,g);ii!=ii_end;++ii) {
            vertex w = source(*ii,g);
            if (v == w) { continue; }
            value_type cache_w = pow(get(wm,*ii),1.0/3.0);
            for (tie(oi,oi_end)=out_edges(w,g);oi!=oi_end;++oi) {
                vertex u = target(*oi,g);
                if (u == w) { continue; }
                if (get(ind,u)) {
                    // cache is cached with its power
                    cur_cc_in += pow(get(wm,*oi),1.0/3.0)*cache_w*cache[u];
                }
            }
        }

        // reset the cache
        for (tie(ii,ii_end)=in_edges(v,g);ii!=ii_end;++ii) {
            if (source(*ii,g) != v) {
                put(ind,source(*ii,g),0);
                cache[source(*ii,g)] = 0;
            }
        }

        // re-init the cache with out-edges
        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            // check to ignore self edges
            if (target(*oi,g) != v) {
                put(ind,target(*oi,g),1);
                cache[target(*oi,g)] += pow(get(wm,*oi),1.0/3.0);
            }
        }

        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            vertex w = target(*oi,g);
            if (v == w) { continue; }
            value_type cache_w = pow(get(wm,*oi),1.0/3.0);
            for (tie(oi2,oi2_end)=out_edges(w,g);oi2!=oi2_end;++oi2) {
                vertex u = target(*oi2,g);
                if (u == w) { continue; }
                if (get(ind,u)) {
                    // cache is cached with its power
                    cur_cc_out += pow(get(wm,*oi2),1.0/3.0)*cache_w*cache[u];
                }
            }
        }

        // store the value
        typename graph_traits<Graph>::degree_size_type
            norm_factor = degs[v]*(degs[v] - 1) - 2*bilateral_edges;
        if (norm_factor > 0) {
            put(cc,v,
                (value_type)(cur_cc_cyc + cur_cc_mid + cur_cc_in + cur_cc_out)/
                    ((value_type)(norm_factor)));
        } else {
            put(cc,v,(value_type)0);
        }

        //std::cout << std::endl;
        //std::cout << v << " " <<  cur_cc_cyc << " " << cur_cc_mid << " " << cur_cc_in << " " << cur_cc_out << " " << degs[v] << " " << bilateral_edges << std::endl;

        for (tie(oi,oi_end)=out_edges(v,g);oi!=oi_end;++oi) {
            // check to ignore self edges
            if (target(*oi,g) != v) {
                put(ind,target(*oi,g),0);
                cache[target(*oi,g)] = 0;
            }
        }
    }
}

int clustering_coefficients(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    double* ccfs, int directed)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_matrix;
    crs_matrix g(nverts, nverts, ia[nverts], ia, ja, NULL);

    std::vector<int> indicator_map(num_vertices(g));

    if (directed) {
        std::vector<mbglIndex> ati(nverts+1);
        std::vector<mbglIndex> atj(ia[nverts]);
        std::vector<mbglIndex> atid(ia[nverts]);

        build_row_and_column_from_csr(g, &ati[0], &atj[0], &atid[0]);

        typedef simple_row_and_column_matrix<mbglIndex,double> bidir_graph;
        bidir_graph bg(nverts, nverts, ia[nverts], ia, ja, NULL, &ati[0], &atj[0], &atid[0]);

        directed_clustering_coefficients(bg,
            make_iterator_property_map(ccfs, get(vertex_index,bg)),
            boost::detail::constant_value_property_map<double>(1.0),
            make_iterator_property_map(indicator_map.begin(),get(vertex_index,bg)));
    } else {
        undirected_clustering_coefficients(g,
            make_iterator_property_map(ccfs, get(vertex_index,g)),
            boost::detail::constant_value_property_map<double>(1.0),
            make_iterator_property_map(indicator_map.begin(), get(vertex_index,g)));
    }

    return (0);
}

int weighted_clustering_coefficients(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double* ccfs)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    std::vector<int> indicator_map(num_vertices(g));

    undirected_clustering_coefficients(g,
        make_iterator_property_map(
            ccfs, get(vertex_index,g)),
        get(edge_weight,g),
        make_iterator_property_map(
            indicator_map.begin(), get(vertex_index,g)));

    return (0);
}

int directed_clustering_coefficients(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double* ccfs)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_matrix;
    crs_matrix g(nverts, nverts, ia[nverts], ia, ja, weight);

    // early return for trivial input
    if (nverts==0) { return (0); }

    std::vector<mbglIndex> ati(nverts+1);
    std::vector<mbglIndex> atj(ia[nverts]);
    std::vector<mbglIndex> atid(ia[nverts]);

    build_row_and_column_from_csr(g, &ati[0], &atj[0], &atid[0]);

    typedef simple_row_and_column_matrix<mbglIndex,double> bidir_graph;
    bidir_graph bg(nverts, nverts, ia[nverts], ia, ja, weight, &ati[0], &atj[0], &atid[0]);

    std::vector<int> indicator_map(num_vertices(g));

    directed_clustering_coefficients(bg,
        make_iterator_property_map(
            ccfs, get(vertex_index,bg)),
        get(edge_weight,bg),
        make_iterator_property_map(
            indicator_map.begin(),
            get(vertex_index,bg)));

    return (0);
}



int betweenness_centrality(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double *centrality, double *ecentrality)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    if (weight)
    {
        if (ecentrality) {
            brandes_betweenness_centrality(g,
	        		centrality_map(centrality)
                    .weight_map(
                        make_iterator_property_map(weight,
                            get(edge_index,g)))
                    .edge_centrality_map(
                        make_iterator_property_map(
                            ecentrality, get(edge_index, g))));
        } else {
            brandes_betweenness_centrality(g,
	        		centrality_map(centrality)
                    .weight_map(
                        make_iterator_property_map(weight,
                            get(edge_index,g))));
        }
    }
    else
    {
        if (ecentrality) {
            brandes_betweenness_centrality(g,
			    	centrality_map(centrality)
                    .edge_centrality_map(
                        make_iterator_property_map(
                            ecentrality, get(edge_index, g))));
        } else {
            brandes_betweenness_centrality(g,
			    	make_iterator_property_map(centrality,
				    	get(vertex_index, g)));
        }
    }

    return (0);
}

/**
 * Test for a topological order or topological sort of a graph.
 *
 * This function will also test if the graph is a dag.  If the
 * rev_order parameter is null, then the actual topological order
 * is ignored and the function just tests for a dag.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param topo_order a topological order of the vertices
 * @param is_dag tests if the graph is a dag (optional)
 */
int topological_order(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex *topo_order, int *is_dag)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    if (is_dag) {
        *is_dag = 1;
    }

    if (topo_order == NULL)
    {
        return (-1);
    }
    else
    {
        typedef reverse_iterator<mbglIndex*>
            reverse_iter;

        // note for the future, boost's reverse iterator is wild in that
        // it actually dereferences the PRIOR element, so the correct
        // range of a boost reverse itertor is the same as the original
        // iterator.  Consequently, the last element is actual topo_order +
        // nverts because the previous element (what * returns) is the
        // correct element.
        reverse_iter output(topo_order + (nverts));

        try {
            topological_sort(g, output);
        } catch (not_a_dag) {
            if (is_dag) {
                *is_dag = 0;
            }
        }
    }

    return (0);
}

template <typename Graph,
        typename MateMap,
        typename VertexIndexMap>
bool matching_help(const Graph& g, MateMap mate, VertexIndexMap vm,
                   int initial_matching, int augmenting_path, int verify)
{
    using namespace boost;
    if (initial_matching == 1)
    {
        if (augmenting_path == 1)
            if (verify == 1)
                return matching<Graph,MateMap,VertexIndexMap,
                        no_augmenting_path_finder,
                        empty_matching,
                        no_matching_verifier>
                        (g, mate, vm);
            else
                return matching<Graph,MateMap,VertexIndexMap,
                        no_augmenting_path_finder,
                        empty_matching,
                        maximum_cardinality_matching_verifier>
                        (g, mate, vm);
        else
            if (verify == 1)
                return matching<Graph,MateMap,VertexIndexMap,
                        edmonds_augmenting_path_finder,
                        empty_matching,
                        no_matching_verifier>
                        (g, mate, vm);
            else
                return matching<Graph,MateMap,VertexIndexMap,
                        edmonds_augmenting_path_finder,
                        empty_matching,
                        maximum_cardinality_matching_verifier>
                        (g, mate, vm);
    }
    else if (initial_matching == 2)
    {
        if (augmenting_path == 1)
            if (verify == 1)
                return matching<Graph,MateMap,VertexIndexMap,
                        no_augmenting_path_finder,
                        greedy_matching,
                        no_matching_verifier>
                        (g, mate, vm);
            else
                return matching<Graph,MateMap,VertexIndexMap,
                        no_augmenting_path_finder,
                        greedy_matching,
                        maximum_cardinality_matching_verifier>
                        (g, mate, vm);
        else
            if (verify == 1)
                return matching<Graph,MateMap,VertexIndexMap,
                        edmonds_augmenting_path_finder,
                        greedy_matching,
                        no_matching_verifier>
                        (g, mate, vm);
            else
                return matching<Graph,MateMap,VertexIndexMap,
                        edmonds_augmenting_path_finder,
                        greedy_matching,
                        maximum_cardinality_matching_verifier>
                        (g, mate, vm);
    }
    else
    {
        if (augmenting_path == 1)
            if (verify == 1)
                return matching<Graph,MateMap,VertexIndexMap,
                        no_augmenting_path_finder,
                        extra_greedy_matching,
                        no_matching_verifier>
                        (g, mate, vm);
            else
                return matching<Graph,MateMap,VertexIndexMap,
                        no_augmenting_path_finder,
                        extra_greedy_matching,
                        maximum_cardinality_matching_verifier>
                        (g, mate, vm);
        else
            if (verify == 1)
                return matching<Graph,MateMap,VertexIndexMap,
                        edmonds_augmenting_path_finder,
                        extra_greedy_matching,
                        no_matching_verifier>
                        (g, mate, vm);
            else
                return matching<Graph,MateMap,VertexIndexMap,
                        edmonds_augmenting_path_finder,
                        extra_greedy_matching,
                        maximum_cardinality_matching_verifier>
                        (g, mate, vm);
    }
}

/** Wrap a boost graph library call to matching.
 *
 * A matching is a subset of edges with no common vertices in an undirected
 * graph.
 *
 * A maximum cardinality matching has the largest number of edges over all
 * possible matchings.
 *
 * The graph must be undirected.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param mate an array of size nverts which stores the matching vertex index
 *  or null vertex if there is no match
 * @param initial_matching indicates the initial matching used in the algorithm
 *  1: no_matching
 *  2: greedy_matching
 *  3: extra_greedy_matching
 * @param augmenting_path indicates the algorithm used to find augmenting paths
 *  1: no_augmenting_path_finder
 *  2: edmonds_augmenting_path_finder
 * @param verify indicates if we verify the algorithm ouput
 *  1: no_matching_verifier
 *  2: maximum_cardinality_matching_verifier
 * @param verified the output of the verification algorithm (if run), or
 *  false otherwise
 * @param null_vertex the special index to indicate an unmatched vertex
 * @return an error code if possible
 *   0: indicates success
 *  -1: indicates a parameter error with initial_matching, augmenting_path, or verify
 */
int maximum_cardinality_matching(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* mate, int initial_matching, int augmenting_path, int verify,
    int *verified, mbglIndex *null_vertex)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    if (initial_matching < 1 || initial_matching > 3 ||
        augmenting_path < 1 || augmenting_path > 2 ||
        verify < 1 || verify > 2)
    {
        return (-1);
    }

    bool max_matching =
        matching_help(g,mate,get(vertex_index,g),
            initial_matching,augmenting_path,verify);

    if (verified) {
        *verified = 0;
        if (verify > 1) { *verified = max_matching; }
    }

    if (null_vertex) {
        *null_vertex = graph_traits<crs_graph>::null_vertex();
    }

    return (0);
}

template <typename Graph, typename MateMap, typename VertexIndexMap>
bool verify_matching_help(const Graph& g, MateMap mate, VertexIndexMap vm)
{
    return boost::maximum_cardinality_matching_verifier<
                    Graph,MateMap,VertexIndexMap>::verify_matching(g,mate,vm);
}

/*
 * mate = nverts if unmatched
 */
int test_maximum_cardinality_matching(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* mate, int *verified)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    // set the null_vertex to the true null
    mbglIndex true_null = graph_traits<crs_graph>::null_vertex();
    for (mbglIndex v=0;v<num_vertices(g);++v) {
        if (mate[v] == nverts) { mate[v] = true_null; }
    }

    bool max_matching = verify_matching_help(g, mate, get(vertex_index,g));
    if (verified) { *verified = 0; }
    if (max_matching && verified) { *verified = 1; }

    return (0);
}

/** Compute the core_numbers of a graph
 *
 * For an undirected graph, this function computes the core number of each
 * vertex.  The core number is the minimum number such that removing all
 * vertices of degree <= cn[k] removes vertex k.  For a directed graph
 * we compute the in-degree core number.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param cn an array of core numbers, length nverts
 * @param rt an array of removal times, length nverts
 * @return an error code if possible
 */
int core_numbers(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglDegreeType *cn, int *rt)
{
    // History
    //
    // 30 July 2007
    // added removal time visitor
    // changed to mbglDegreeType
    using namespace yasmic;
    using namespace boost;

    if (!cn || !rt) { return (-1); }


    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    int time = 0;

    core_numbers(g, make_iterator_property_map(cn, get(vertex_index,g)),
        make_core_numbers_visitor(stamp_times(rt, time, on_examine_vertex())));

    return (0);
}

/** Compute the weighted core numbers of a graph
 *
 * For an undirected graph, this function computes the core number of each
 * vertex.  The core number is the minimum number such that removing all
 * vertices of weighted in-degree <= cn[k] removes vertex k.  For a
 * directed graph we compute the weighted in-degree core number.
 *
 * @param nverts the number of vertices in the graph
 * @param ja the connectivity for each vertex
 * @param ia the row connectivity points into ja
 * @param cn an array of core numbers, length nverts
 * @param rt an array of removal times, length nverts
 * @return an error code if possible
 */
int weighted_core_numbers(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double *cn, int *rt)
{
    // History
    //
    // 30 July 2007
    // added removal time visitor
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    int time=0;

    weighted_core_numbers(g, make_iterator_property_map(cn, get(vertex_index,g)),
        make_core_numbers_visitor(stamp_times(rt, time, on_examine_vertex())));

    return (0);
}


int dominator_tree(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    mbglIndex src, mbglIndex *pred)
{
    //
    // History
    //
    // 12 July 2007
    // Initial implementation
    //
    // 23 July 2007
    // Added iterator property map call for pred.
    //
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_graph;
    crs_graph g(nverts, nverts, ia[nverts], ia, ja, NULL);

    std::vector<mbglIndex> ati(nverts+1);
    std::vector<mbglIndex> atj(ia[nverts]);
    std::vector<mbglIndex> atid(ia[nverts]);

    build_row_and_column_from_csr(g, &ati[0], &atj[0], &atid[0]);

    typedef simple_row_and_column_matrix<mbglIndex,double> bidir_graph;
    bidir_graph bg(nverts, nverts, ia[nverts], ia, ja, NULL, &ati[0], &atj[0], &atid[0]);

    // modify the output to conform to what matlab_bgl expects
    mbglIndex null_vertex = graph_traits<bidir_graph>::null_vertex();
    for (mbglIndex i=0; i< nverts; i++) { pred[i] = null_vertex; }

    lengauer_tarjan_dominator_tree(bg, src,
        make_iterator_property_map(pred, get(vertex_index,bg)));

    pred[src] = src;

    // modify the output to conform to what matlab_bgl expects
    for (mbglIndex i=0; i< nverts; i++) {
        if (pred[i] == null_vertex) { pred[i] = i; }
    }

    return (0);
}

