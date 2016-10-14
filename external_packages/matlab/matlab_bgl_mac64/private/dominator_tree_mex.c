/*
 * David Gleich
 * Copyright, Stanford Unviersity, 2007
 */

/** 
 * @file dominator_tree_mex.c
 * Wrap a call to the libmbgl shortest path functions.
 */

/*
 * 12 July 2007
 * Initial version
 */


#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"
#include "expand_macros.h"
#include "common_functions.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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
    
    /* source */
    mwIndex u;
   
    /* output data */
    double *pred;
    
    /* 
     * The current calling pattern is
     * dominator_tree_mex(A,u)
     * u is the source vertex
     */
    
    const mxArray* arg_matrix;
    const mxArray* arg_source;
    int required_arguments = 2;
    
    if (nrhs != required_arguments) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the function requires %i arguments, not %i\n", 
            required_arguments, nrhs);
    }
    
    arg_matrix = prhs[0];
    arg_source = prhs[1];
    
    u = (mwIndex)load_scalar_arg(arg_source,1);
    
    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    if (mrows != ncols || !mxIsSparse(arg_matrix)) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the matrix must be sparse and square");
    }

    n = mrows;
    
    /* Get the sparse matrix */
    
    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);
    
    nz = ia[n];
    
    /* Process the scalar target */
    u = u-1;
    
    if (u < 0 || u >= n) 
    {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexParameter", 
            "start vertex (%i) not a valid vertex.", u+1);
    }
    
    
    plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
    
    /* create the output vectors */
    
    pred = mxGetPr(plhs[0]);
    
    
    #ifdef _DEBUG
    mexPrintf("dominator_tree...");
    #endif 
    
    dominator_tree(n, ja, ia, u, (mwIndex*)pred);
    
    #ifdef _DEBUG
    mexPrintf("done\n");
    #endif 
    
    expand_index_to_double_zero_equality((mwIndex*)pred, pred, n, 1.0);
    /* expand_index_to_double((mwIndex*)pred, pred, n, 0.0); */
   
    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif 
}


