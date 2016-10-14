/*
 * ==============================================================
 * bfs_mex.c The mex interface to the matlab bgl wrapper.
 *
 * David Gleich
 * 20 April 20020
 * =============================================================
 */

/*
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 *
 * 1 March 2007
 * Updated to use expand macros
 *
 * 19 April 2007
 * Fixed error with invalid vertex by checking the input to make sure the
 * start vertex is valid.
 *
 * Added target vertex
 */


#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"

#include <math.h>
#include <stdlib.h>

#include "expand_macros.h"

/*
 * The mex function runs a max-flow min-cut problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwIndex i;
    
    mwIndex mrows, ncols;
    
    mwIndex n,nz;
    
    /* sparse matrix */
    mwIndex *ia, *ja;
    
    /* source */
    mwIndex u;
    
    /* target */
    mwIndex v;
    int full;
    
    /* output data */
    double *d, *dt, *ft, *pred;
    
    int *int_d;
    int *int_dt;
    int *int_ft;
    
    /* 
     * The current calling pattern is
     * dfs_mex(A,u,v,full)
     * 
     * u and v are the source and target.  v = 0 if there is no target 
     * and the entire search should complete.
     */
    
    const mxArray* arg_matrix;
    const mxArray* arg_source;
    const mxArray* arg_target;
    const mxArray* arg_full;
    
    if (nrhs != 4) 
    {
        mexErrMsgTxt("4 inputs required.");
    }
    
    arg_matrix = prhs[0];
    arg_source = prhs[1];
    arg_target = prhs[2];
    arg_full = prhs[3];

    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    if (mrows != ncols ||
        !mxIsSparse(arg_matrix)) 
    {
        mexErrMsgTxt("Input must be a noncomplex square sparse matrix.");
    }
    
    n = mrows;
        
    /* The second input must be a scalar. */
    if (mxGetNumberOfElements(arg_source) > 1 || !mxIsDouble(arg_source))
    {
        mexErrMsgTxt("Invalid scalar.");
    }
    
    /* Get the sparse matrix */
    
    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);
    
    nz = ia[n];
    
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
    
    full = (int)mxGetScalar(arg_full);
    
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,n,mxREAL);
    
    d = mxGetPr(plhs[0]);
    dt = mxGetPr(plhs[1]);
    ft = mxGetPr(plhs[2]);
    pred = mxGetPr(plhs[3]);
    
    int_d = (int*)d;
    int_dt = (int*)dt;
    int_ft = (int*)ft;
    
    for (i=0; i < n; i++)
    {
        int_d[i]=-1;
        int_dt[i]=-1;
        int_ft[i]=-1;
    }
    
    int_d[u] = 0;
    
    #ifdef _DEBUG
    mexPrintf("dfs...");
    #endif 
    depth_first_search(n, ja, ia,
        u, full, v, (int*)d, (int*)dt, (int*)ft, (mwIndex*)pred);
    #ifdef _DEBUG
    mexPrintf("done!\n");
    #endif 
   
    
    expand_int_to_double((int*)d, d, n, 0.0);
    expand_int_to_double((int*)dt, dt, n, 0.0);
    expand_int_to_double((int*)ft, ft, n, 0.0);
    expand_index_to_double_zero_equality((mwIndex*)pred, pred, n, 1.0);
    
    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif 
}

