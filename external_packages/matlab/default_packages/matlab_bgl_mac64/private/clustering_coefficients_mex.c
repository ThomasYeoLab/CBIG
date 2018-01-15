/** @file clustering_coefficients_mex.c
 * @copyright Stanford University, 2006-2008
 * @author David F. Gleich
 * The mex interface to the matlab bgl wrapper.
 */

/** History
 *  2006-04-23: Initial version
 *  2007-02-19: Updated to Matlab 2006b sparse matrix interface
 *  2007-07-11: Update for weighted and directed clusternig coefficients
 */

#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"

#include "common_functions.h"

#include <math.h>
#include <stdlib.h>

/*
 * The mex function runs a clustering coefficients problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwIndex mrows, ncols;

    mwIndex n,nz;

    /* sparse matrix */
    mwIndex *ia, *ja;
    double *a;

    /* output data */
    double *ccfs;

    /* used to switch between algorithm types */
    int weight_type; /* = 0 if unweighted,
                        = 1 if a vector of weights,
                        = 2 if given by the matrix */
    int undirected;

    /*
     * The current calling pattern is
     * clustering_coefficients_mex(A,undirected,weight)
     * where weight = 0 to use the unweighted version
     *       weight = 'matrix' to use the values in the matrix
     *       weight = vector to use a vector of weights
     */

    const mxArray* arg_matrix;
    const mxArray* arg_undirected;
    const mxArray* arg_weight;
    int required_arguments = 3;

    if (nrhs != required_arguments) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the function requires %i arguments, not %i\n",
            required_arguments, nrhs);
    }

    arg_matrix = prhs[0];
    arg_undirected = prhs[1];
    arg_weight = prhs[2];

    undirected = (int)load_scalar_arg(arg_undirected, 1);

    if (mxGetNumberOfElements(arg_weight) == 1)
    {
        /* make sure it is valid */
        if ((int)mxGetScalar(arg_weight) != 0) {
            mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
                "unknown weight option %g\n", mxGetScalar(arg_weight));
        }

        weight_type = 0;
        a = NULL;
    }
    else if (mxIsChar(arg_weight)) {
        weight_type = 2;
        a = mxGetPr(arg_matrix);
    }
    else if (mxIsDouble(arg_weight)) {
        weight_type = 1;
        a = mxGetPr(arg_weight);
    }
    else {
        mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
            "unrecognized weight option");
        return;
    }

    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    if (mrows != ncols || !mxIsSparse(arg_matrix) ||
        ((!mxIsDouble(arg_matrix) || mxIsComplex(arg_matrix)) && (weight_type == 2))
        )
    {
        mexErrMsgTxt("Input must be a square sparse matrix.");
    }

    n = mrows;

    /* recall that we've transposed the matrix */
    ja = mxGetIr(prhs[0]);
    ia = mxGetJc(prhs[0]);

    nz = ia[n];

    /* check the reweighting array to make sure it is acceptable */
    if (weight_type == 1 && mxGetNumberOfElements(arg_weight) < nz) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
            "the weight array must have length >= nnz(A)");
    }

    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);

    ccfs = mxGetPr(plhs[0]);

    #ifdef _DEBUG
    mexPrintf("clustering_coefficients...");
    #endif

    if (weight_type == 0) {
        clustering_coefficients(n, ja, ia, ccfs, !undirected);
    } else if (undirected) {
        weighted_clustering_coefficients(n, ja, ia, a, ccfs);
    } else {
        directed_clustering_coefficients(n, ja, ia, a, ccfs);
    }


    #ifdef _DEBUG
    mexPrintf("done\n");
    #endif

    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}

