/*
 * David Gleich
 * Copyright, Stanford Unviersity, 2007
 */

/**
 * @file matlab_bgl_sp_mex.c
 * Wrap a call to the libmbgl shortest path functions.
 */

/*
 * 20 April 2006
 * Initial version
 *
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 *
 * 1 March 2007
 * Updated to use expand macros
 *
 * 18 April 2007
 * Updated to support additional 'weight' parameter, this involved
 * changing the calling pattern of the function such that the visitor
 * would be the 6th option and not the 5th.
 *
 * 19 April 2007
 * Fixed error with invalid vertex by checking the input to make sure the
 * start vertex is valid.
 *
 * 12 July 2007
 * Updated header information
 * Updated to use load_string_arg
 */


#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"
#include "visitor_macros.h"
#include "expand_macros.h"
#include "common_functions.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

PROTOTYPE_VISITOR_VERTEX_FUNCTION(initialize_vertex)
PROTOTYPE_VISITOR_VERTEX_FUNCTION(examine_vertex)
PROTOTYPE_VISITOR_VERTEX_FUNCTION(discover_vertex)
PROTOTYPE_VISITOR_VERTEX_FUNCTION(finish_vertex)

PROTOTYPE_VISITOR_EDGE_FUNCTION(examine_edge)
PROTOTYPE_VISITOR_EDGE_FUNCTION(edge_relaxed)
PROTOTYPE_VISITOR_EDGE_FUNCTION(edge_not_relaxed)
PROTOTYPE_VISITOR_EDGE_FUNCTION(edge_minimized)
PROTOTYPE_VISITOR_EDGE_FUNCTION(edge_not_minimized)

/*
 * The mex function runs a shortest path problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwIndex mrows, ncols;

    mwIndex n,nz;

    /* sparse matrix */
    mwIndex *ia, *ja;
    double *a;

    /* source */
    mwIndex u;

    /* target */
    mwIndex v;

    char *algname;

    double dinf;

    /* true if this function was called with a visitor */
    int use_visitor = 0;

    /* true if this function is reweighted */
    int reweighted = 0;

    /* output data */
    double *d, *pred;

    /*
     * The current calling pattern is
     * matlab_bgl_sp_mex(A,u,v,algname,dinf,reweight,[visitor])
     * so visitor is an optional paramter.
     * algname is a string with either 'dag', 'dijkstra' or 'bellman_ford'
     * reweight is either a string of a length nnz vector
     *
     * if reweight is a length nnz vector, then we use that as the values
     * for the matrix, if its a string, then we use the values from the
     * the matrix.
     *
     * u and v are the source and target.  v = 0 if there is no target
     * and the entire search should complete.
     */

    const mxArray* arg_matrix;
    const mxArray* arg_source;
    const mxArray* arg_target;
    const mxArray* arg_algname;
    const mxArray* arg_dinf;
    const mxArray* arg_reweight;
    const mxArray* arg_visitor=NULL;

    if (nrhs < 6 || nrhs > 7)
    {
        mexErrMsgTxt("6 or 7 inputs required.");
    }

    arg_matrix = prhs[0];
    arg_source = prhs[1];
    arg_target = prhs[2];
    arg_algname = prhs[3];
    arg_dinf = prhs[4];
    arg_reweight = prhs[5];

    algname = load_string_arg(arg_algname,3);

    /* First test if they are going to reweight or not */
    if (!mxIsChar(arg_reweight))
    {
        reweighted = 1;
    }


    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);

    if (!reweighted && (mrows != ncols ||
        !mxIsSparse(arg_matrix) ||
        !mxIsDouble(arg_matrix) ||
        mxIsComplex(arg_matrix)))
    {
        mexErrMsgTxt("A non-reweighted input matrix must be a noncomplex square sparse matrix.");
    }
    else if (reweighted && (mrows != ncols ||
        !mxIsSparse(arg_matrix)))
    {
        mexErrMsgTxt("A reweighted input matrix must be a square sparse matrix.");
    }

    /* The 6th input (if present) must be a structure. */
    if (nrhs == 7 && !mxIsStruct(prhs[6]))
    {
        mexErrMsgTxt("Invalid structure.");
    }

    if (nrhs == 7)
    {
        arg_visitor = prhs[6];
        use_visitor = 1;
    }

    n = mrows;


    /* Get the sparse matrix */

    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);

    nz = ia[n];

    if (!reweighted)
    {
        a = mxGetPr(arg_matrix);
    }
    else
    {
        /* check the reweighting array to make sure it is acceptable */
        if (mxGetNumberOfElements(arg_reweight) < nz)
            mexErrMsgTxt("The reweight array must have length at least nnz(A)");
        if (!mxIsDouble(arg_reweight))
            mexErrMsgTxt("The reweight array must be of type double.");

        a = mxGetPr(arg_reweight);
    }

    /* Get the scalar */
    u = (mwIndex)mxGetScalar(arg_source);
    u = u-1;

    /* Get the target */
    v = (mwIndex)mxGetScalar(arg_target);
    if (v == 0) {
        /* they want the entire search */
        v = n;
    } else {
        /* they want to look for vertex v */
        v -= 1;
    }

    if (u < 0 || u >= n)
    {
        mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
            "start vertex (%i) not a valid vertex.", u+1);
    }

    if (v < 0 || v > n)
    {
        mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
            "target vertex (%i) not a valid vertex.", v+1);
    }


    /* Get the uninitialized value */
    dinf = mxGetScalar(arg_dinf);

    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,n,mxREAL);

    /* create the output vectors */

    d = mxGetPr(plhs[0]);
    pred = mxGetPr(plhs[1]);


    #ifdef _DEBUG
    mexPrintf("sp_%s...",algname);
    #endif
    if (strcmp(algname, "dijkstra") == 0)
    {
        if (use_visitor)
        {
            const mxArray *vis = arg_visitor;
            dijkstra_visitor_funcs_t d_vis = {0};

            /* Check the visitor and construct the visitor structure. */
            d_vis.pdata = (void*)vis;
            CHECK_AND_SET_VISITOR_FUNCTION(vis,initialize_vertex,d_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,discover_vertex,d_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_vertex,d_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,finish_vertex,d_vis);

            CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_edge,d_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_relaxed,d_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_not_relaxed,d_vis);

            dijkstra_sp_visitor(n, ja, ia, a,
                u,
                d, (mwIndex*)pred, dinf, d_vis);
        }
        else
        {
            dijkstra_sp(n, ja, ia, a,
                u, v,
                d, (mwIndex*)pred, dinf);
        }
    }
    else if (strcmp(algname, "bellman_ford") == 0)
    {
        if (use_visitor)
        {
            const mxArray *vis = arg_visitor;
            bellman_ford_visitor_funcs_t bf_vis = {0};

            /* Check the visitor and construct the visitor structure. */
            bf_vis.pdata = (void*)vis;
            CHECK_AND_SET_VISITOR_FUNCTION(vis,initialize_vertex,bf_vis);

            CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_edge,bf_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_relaxed,bf_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_not_relaxed,bf_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_minimized,bf_vis);
            CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_not_minimized,bf_vis);

            bellman_ford_sp_visitor(n, ja, ia, a,
            u,
            d, (mwIndex*)pred, dinf, bf_vis);
        }
        else
        {
            bellman_ford_sp(n, ja, ia, a,
            u, v,
            d, (mwIndex*)pred, dinf);
        }

    }
    else if (strcmp(algname, "dag") == 0)
    {
        if (use_visitor) { mexWarnMsgTxt("Visitor ignored."); }
        dag_sp(n, ja, ia, a,
            u, v,
            d, (mwIndex*)pred, dinf);
    }
    else
    {
        mexErrMsgTxt("Unknown algname.");
    }
    #ifdef _DEBUG
    mexPrintf("done\n");
    #endif

    expand_index_to_double_zero_equality((mwIndex*)pred, pred, n, 1.0);

    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}



CALL_MATLAB_VERTEX_VISITOR_FUNCTION(initialize_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(examine_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(discover_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(finish_vertex)

CALL_MATLAB_EDGE_VISITOR_FUNCTION(examine_edge)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(edge_relaxed)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(edge_not_relaxed)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(edge_minimized)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(edge_not_minimized)

