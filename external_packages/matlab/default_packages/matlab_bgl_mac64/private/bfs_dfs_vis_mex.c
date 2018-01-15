/** @file bfs_dfs_vis_mex.c
 * @author David F. Gleich
 * @date 2008-09-29
 * @copyright Stanford University, 2006-2008
 * The mex interface to the matlab bgl wrapper for bfs and dfs with a visitor.
 */

/** History
 *  2006-05-30: Initial version
 *  2007-02-19: Updated to use Matlab 2006b sparse matrix interface
 *  2007-02-22: Updated to use large graph libmbgl array
 *  2007-04-19: Fixed error with invalid vertex by checking the
 *    input to make sure the start vertex is valid.
 */


#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"
#include "visitor_macros.h"
#include <math.h>
#include <stdlib.h>

int call_matlab_initialize_vertex(void *pdata, mwIndex u);
int call_matlab_discover_vertex(void *pdata, mwIndex u);
int call_matlab_examine_vertex(void *pdata, mwIndex u);
int call_matlab_finish_vertex(void *pdata, mwIndex u);
int call_matlab_start_vertex(void *pdata, mwIndex u);

int call_matlab_examine_edge(void *pdata, mwIndex ei, mwIndex u, mwIndex v);
int call_matlab_tree_edge(void *pdata, mwIndex ei, mwIndex u, mwIndex v);
int call_matlab_non_tree_edge(void *pdata, mwIndex ei, mwIndex u, mwIndex v);
int call_matlab_gray_target(void *pdata, mwIndex ei, mwIndex u, mwIndex v);
int call_matlab_black_target(void *pdata, mwIndex ei, mwIndex u, mwIndex v);
int call_matlab_back_edge(void *pdata, mwIndex ei, mwIndex u, mwIndex v);
int call_matlab_forward_or_cross_edge(void *pdata, mwIndex ei, mwIndex u, mwIndex v);

void bfs_vis(mwIndex n, mwIndex *ja, mwIndex *ia, mwIndex u,
    const mxArray *vis)
{
    bfs_visitor_funcs_t bfs_vis = {0};

    /* Check the visitor and construct the visitor structure. */
    bfs_vis.pdata = (void*)vis;

    CHECK_AND_SET_VISITOR_FUNCTION(vis,initialize_vertex,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,discover_vertex,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_vertex,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,finish_vertex,bfs_vis);

    CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_edge,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,tree_edge,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,non_tree_edge,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,gray_target,bfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,black_target,bfs_vis);

    #ifdef _DEBUG
    mexPrintf("bfs...");
    #endif

    breadth_first_search_visitor(n, ja, ia, u, bfs_vis);
}

void dfs_vis(mwIndex n, mwIndex *ja, mwIndex *ia, mwIndex u,
    int full, const mxArray *vis)
{
    dfs_visitor_funcs_t dfs_vis = {0};

    /* Check the visitor and construct the visitor structure. */
    dfs_vis.pdata = (void*)vis;
    CHECK_AND_SET_VISITOR_FUNCTION(vis,initialize_vertex,dfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,discover_vertex,dfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,start_vertex,dfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,finish_vertex,dfs_vis);

    CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_edge,dfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,tree_edge,dfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,back_edge,dfs_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,forward_or_cross_edge,dfs_vis);

    #ifdef _DEBUG
    mexPrintf("dfs...");
    #endif

    depth_first_search_visitor(n, ja, ia, u, full, dfs_vis);
}

/*
 * The mex function runs a BFS or DFS with a visitor.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwIndex mrows, ncols;
    mwIndex n,nz;

    /* sparse matrix */
    mwIndex *ia, *ja;

    /* start */
    mwIndex u;

    const mxArray *vis;

    int call;

    if (nrhs != 4)
    {
        mexErrMsgTxt("2 inputs required.");
    }

    /* The first input must be a sparse matrix. */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if (mrows != ncols ||
        !mxIsSparse(prhs[0]))
    {
        mexErrMsgTxt("Input must be a square sparse matrix.");
    }

    n = mrows;

    /* The second input must be a scalar. */
    if (mxGetNumberOfElements(prhs[1]) > 1 || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("Invalid scalar.");
    }

    /* The third input must be a structure. */
    if (!mxIsStruct(prhs[2]))
    {
        mexErrMsgTxt("Invalid structure.");
    }


    /* The fourth input must be a scalar. */
    if (mxGetNumberOfElements(prhs[3]) > 1 || !mxIsDouble(prhs[3]))
    {
        mexErrMsgTxt("Invalid scalar.");
    }


    /* Get the sparse matrix */

    /* recall that we've transposed the matrix */
    ja = mxGetIr(prhs[0]);
    ia = mxGetJc(prhs[0]);

    nz = ia[n];

    /* Get the scalar */
    u = (mwIndex)mxGetScalar(prhs[1]);
    u = u-1;

    if (u < 0 || u >= n)
    {
        mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
            "start vertex (%i) not a valid vertex.", u+1);
    }

    vis = prhs[2];

    call = (int)mxGetScalar(prhs[3]);

    if (call == 101)
    {
        bfs_vis(n, ja, ia, u, vis);
    }
    else if (call == 201)
    {
        /* call without the "full" option */
        dfs_vis(n, ja, ia, u, 0, vis);
    }
    else if (call == 202)
    {
        /* call with the "full" option */
        dfs_vis(n, ja, ia, u, 1, vis);
    }
    else
    {
        mexErrMsgTxt("Invalid call.");
    }
    #ifdef _DEBUG
    mexPrintf("done!\n");
    #endif


    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}


CALL_MATLAB_VERTEX_VISITOR_FUNCTION(initialize_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(discover_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(examine_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(finish_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(start_vertex)

CALL_MATLAB_EDGE_VISITOR_FUNCTION(examine_edge)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(tree_edge)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(non_tree_edge)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(gray_target)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(black_target)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(back_edge)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(forward_or_cross_edge)

