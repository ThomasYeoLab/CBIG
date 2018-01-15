/** @file core_numbers_mex.c
 * @copyright Stanford University, 2007-2008
 * @author David F. Gleich
 * Wrap a call to the libmbgl core_numbers function.
 */

/** History
 *  2007-07-08: Initial version
 *  2007-07-11: Updated for weighted cores
 *  2007-07-30: Added option for removal times
 */

#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"

#include "expand_macros.h"

#include <math.h>
#include <stdlib.h>

/*
 * The mex function runs a connected components problem.
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
    double *cn;
    double *rt;


    /* used to switch between algorithm types */
    int weight_type; /* = 0 if unweighted,
                        = 1 if a vector of weights,
                        = 2 if given by the matrix */

    /*
     * The current calling pattern is
     * core_numbers_mex(A,weight)
     * where weight = 0 to use the unweighted version
     *       weight = 'matrix' to use the values in the matrix
     *       weight = vector to use a vector of weights
     */

    const mxArray* arg_matrix;
    const mxArray* arg_weight;
    int required_arguments = 2;

    if (nrhs != required_arguments) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the function requires %i arguments, not %i\n",
            required_arguments, nrhs);
    }

    arg_matrix = prhs[0];
    arg_weight = prhs[1];

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
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);

    cn = mxGetPr(plhs[0]);
    rt = mxGetPr(plhs[1]);

    #ifdef _DEBUG
    mexPrintf("core_numbers...");
    #endif

    if (weight_type == 0) {
        core_numbers(n, ja, ia, (mwIndex*)cn, (int*)rt);
    } else {
        weighted_core_numbers(n, ja, ia, a, cn, (int*)rt);
    }


    #ifdef _DEBUG
    mexPrintf("done!\n");
    #endif

    if (weight_type == 0) {
        expand_index_to_double((mwIndex*)cn,cn,n,0.0);
    }
    expand_int_to_double((int*)rt,rt,n,0.0);

    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}



