/** @file betweenness_centrality_mex.c
 * @copyright Stanford University, 2007-2008
 * @author David F. Gleich
 * The mex interface to the matlab bgl wrapper for the topological order func
 */

/** History
 *  2007-04-20: Initial version
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
 * The mex function runs a topological order problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int rval;
    int mrows, ncols;

    int n,nz;

    /* sparse matrix */
    mwIndex *ia, *ja;

    /* output data */
    double *perm;
    double *dag;

    /* arguments */
    const mxArray* arg_matrix;

    if (nrhs != 1)
    {
        mexErrMsgTxt("1 inputs required.");
    }

    arg_matrix = prhs[0];

    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    if (mrows != ncols ||
        !mxIsSparse(arg_matrix))
    {
        mexErrMsgTxt("Input must be a noncomplex square sparse matrix.");
    }

    n = mrows;


    /* Get the sparse matrix */

    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);

    nz = ia[n];

    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

    perm = mxGetPr(plhs[0]);
    dag = mxGetPr(plhs[1]);

    #ifdef _DEBUG
    mexPrintf("topological_order...");
    #endif

    rval = topological_order(n, ja, ia, (mwIndex*)perm, (int*)dag);

    #ifdef _DEBUG
    mexPrintf("done, rval=%i\n", rval);
    #endif

    expand_index_to_double((mwIndex*)perm, perm, n, 1.0);

    expand_int_to_double((int*)dag, dag, 1, 0.0);

    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}

