/*
 * David F. Gleich
 * Copyright, Stanford University, 2006-2008
 */

/** History
 *  2007-07-04: Added maximum_cardinality_matching prototype
 *    Added dominator_tree prototype
 *  2008-09-25: Added kamada_kawai_spring_layout prototype
 *    Added fruchterman_reingold_force_directed_layout prototype
 *    Added gursoy_atun_layout prototype
 *  2008-10-06: Removed prototype comments from this file
 */

#ifndef MATLAB_BGL_H
#define MATLAB_BGL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "matlab_bgl_types.h"

#define MATLAB_BGL_VISITOR_METHOD_VERTEX(a) int (*a)(void *pdata, mbglIndex u)
#define MATLAB_BGL_VISITOR_METHOD_EDGE(a) int (*a)(void *pdata, mbglIndex ei, mbglIndex u, mbglIndex v)

typedef struct {
    void *pdata;

    MATLAB_BGL_VISITOR_METHOD_VERTEX(initialize_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(discover_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(examine_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(finish_vertex);

    MATLAB_BGL_VISITOR_METHOD_EDGE(examine_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(tree_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(non_tree_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(gray_target);
    MATLAB_BGL_VISITOR_METHOD_EDGE(black_target);

} bfs_visitor_funcs_t;

typedef struct {
    void *pdata;

    MATLAB_BGL_VISITOR_METHOD_VERTEX(initialize_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(start_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(discover_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(finish_vertex);

    MATLAB_BGL_VISITOR_METHOD_EDGE(examine_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(tree_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(back_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(forward_or_cross_edge);
} dfs_visitor_funcs_t;

typedef struct {
    void *pdata;

    MATLAB_BGL_VISITOR_METHOD_VERTEX(initialize_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(examine_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(discover_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(finish_vertex);

    MATLAB_BGL_VISITOR_METHOD_EDGE(examine_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_relaxed);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_not_relaxed);
} dijkstra_visitor_funcs_t;

typedef struct {
    void *pdata;

    MATLAB_BGL_VISITOR_METHOD_VERTEX(initialize_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(examine_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(discover_vertex);
    MATLAB_BGL_VISITOR_METHOD_VERTEX(finish_vertex);

    MATLAB_BGL_VISITOR_METHOD_EDGE(examine_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_relaxed);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_not_relaxed);
    MATLAB_BGL_VISITOR_METHOD_EDGE(black_target);
} astar_visitor_funcs_t;


typedef struct {
    void *pdata;

    MATLAB_BGL_VISITOR_METHOD_VERTEX(initialize_vertex);

    MATLAB_BGL_VISITOR_METHOD_EDGE(examine_edge);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_relaxed);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_not_relaxed);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_minimized);
    MATLAB_BGL_VISITOR_METHOD_EDGE(edge_not_minimized);

} bellman_ford_visitor_funcs_t;

/**
 * @section max_flow.cc
 * Prototypes for max_flow functions
 */

int push_relabel_max_flow(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,  /* connectivity params */
    mbglIndex src, mbglIndex sink, /* flow data */
    int* cap, int* res, /* capacity and residual capacity */
    mbglIndex* rev_edge_index, int *flow);

int kolmogorov_max_flow(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,  /* connectivity params */
    mbglIndex src, mbglIndex sink, /* flow data */
    int* cap, int* res, /* capacity and residual capacity */
    mbglIndex* rev_edge_index, int *flow);

int edmunds_karp_max_flow(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,  /* connectivity params */
    mbglIndex src, mbglIndex sink, /* flow data */
    int* cap, int* res, /* capacity and residual capacity */
    mbglIndex* rev_edge_index, int *flow);

/**
 * @section searches.cc
 */

int breadth_first_search(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    int* d, int* dt, mbglIndex* pred /* output data: distance, discover time, predecessor */
    );

int breadth_first_search_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, /* problem data */
    bfs_visitor_funcs_t vis /* visitor structure */
    );

int depth_first_search(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, int full, mbglIndex dst, /* problem data */
    int* d, int* dt, int* ft, mbglIndex *pred /* output data: discover time, finish time, predecessor */
    );

int depth_first_search_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex src, int full, /* problem data */
    dfs_visitor_funcs_t vis /* visitor structure */
    );

int astar_search(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* start nd target vertices */
    double *d, mbglIndex *pred, double *f, /* output */
    double *h /* heuristic function value for all vertices */, double dinf);

int astar_search_hfunc(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* start and target vertices */
    double *d, mbglIndex *pred, double *f,
    double (*hfunc)(void* pdata, mbglIndex u), void* pdata /* heuristic function */, double dinf);

int astar_search_hfunc_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, /* start vertex */
    double *d, mbglIndex *pred, double *f,
    double (*hfunc)(void* pdata, mbglIndex u), void* pdata /* heuristic function */, double dinf,
    astar_visitor_funcs_t vis);

/**
 * @section components.cc
 */

int biconnected_components(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* a, mbglIndex* ci);

int strong_components(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* ci);

/**
 * @section shortest_path.cc
 */

int dijkstra_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf);

int dijkstra_sp_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, /* problem data */
    double* d, mbglIndex *pred, double dinf,
    dijkstra_visitor_funcs_t vis);

int bellman_ford_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf);

int bellman_ford_sp_visitor(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, /* problem data */
    double* d, mbglIndex *pred, double dinf,
    bellman_ford_visitor_funcs_t vis);

int dag_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf);

int johnson_all_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double* D, double dinf);

int floyd_warshall_all_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double* D, double dinf, mbglIndex *pred);

/**
 * @section spanning_trees.cc
 */

int kruskal_mst(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex* i, mbglIndex* j, double* val, mbglIndex* nedges /* tree output */);

int prim_mst(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex* i, mbglIndex* j, double* val, mbglIndex* nedges /* tree output */);

int prim_mst_rooted(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex* i, mbglIndex* j, double* val, mbglIndex* nedges, /* tree output */
    mbglIndex root);

/**
 * @section statistics.cc
 *
 */
int betweenness_centrality(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double *centrality, double *ecentrality);

int clustering_coefficients(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    double *ccoeffs, int directed);

int weighted_clustering_coefficients(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double *ccoeffs);

int directed_clustering_coefficients(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double *ccoeffs);

int topological_order(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex *rev_order, int *is_dag);

int maximum_cardinality_matching(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* mate, int initial_matching, int augmenting_path, int verify,
    int *verified, mbglIndex *null_vertex);

int test_maximum_cardinality_matching(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex* mate, int *verified);

int core_numbers(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglDegreeType *cn, int *rt);

int weighted_core_numbers(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    double *cn, int *rt);

int dominator_tree(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    mbglIndex src, mbglIndex *pred);

/**
 * @section layouts.cc
 * Prototypes for layouts.cc
 */

int kamada_kawai_spring_layout(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight,
    double tol, int maxiter, double spring_constant, int progressive, 
    double edge_length,
    double *positions,
    double *spring_strength, double *distance);

int fruchterman_reingold_force_directed_layout(
		mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
		int iterations, double initial_temp, int grid_force_pairs,
		double width, double height, int progressive,
		double *positions);

/** The possible topologies for the gursoy_atun_layout
 */
typedef enum  {
	CUBE_LAYOUT_TOPOLOGY, BALL_LAYOUT_TOPOLOGY, HEART_LAYOUT_TOPOLOGY
} gursoy_atun_layout_topology_t;

extern const int gursoy_atun_invalid_dim;
extern const int gursoy_atun_dim_too_large;

int gursoy_atun_layout(
		mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight,
		gursoy_atun_layout_topology_t topology, int topology_dim,
		int iterations, double diameter_i, double diameter_f,
		double learning_constant_i, double learning_constant_f,
		double *positions);

/**
 * @section planar.cc
 * Prototypes for filse in planar.cc
 */
int boyer_myrvold_planarity_test(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int *is_planar,
    mbglIndex *i, mbglIndex *j, mbglIndex* nedges, /* kuratowski subgraph output */
    mbglIndex *eip, mbglIndex *eie); /* embedding information */

int chrobak_payne_straight_line_drawing(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int just_ordering, int is_maximal,
    mbglIndex *i, mbglIndex *j, mbglIndex* nedges, /* extra edges */
    mbglIndex *p, /* ordering permutation */
    mbglDegreeType *X);

int is_kuratowski_subgraph(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int *is_ksubgraph);

int is_straight_line_drawing(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    double *X, int *is_sldrawing);

int triangulate_graph(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia,
    int make_connected, int make_biconnected, int make_maximal,
    mbglIndex *i, mbglIndex *j, mbglIndex* nedges /* extra edges */);

/**
 * @section orderings.cc
 */

int reverse_cuthill_mckee_order(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex start, /* input */
    mbglIndex *perm /* permutation output */);

int king_order(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex *perm /* permutation output */);

int minimum_degree_order(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex *perm /* permutation output */);

int sloan_order(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, /* connectivity params */
    mbglIndex *perm /* permutation output */);

#ifdef __cplusplus
}
#endif /* __cplusplus */



#endif /* MATLAB_BGL_H */


