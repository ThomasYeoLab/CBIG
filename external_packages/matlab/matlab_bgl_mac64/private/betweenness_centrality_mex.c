/** @file betweenness_centrality_mex.c
 * @copyright Stanford University, 2006-2008
 * @author David F. Gleich
 * The mex interface to the matlab bgl wrapper.
 */

/** History
 *  2006-04-23: Initial version
 *  2007-02-19: Updated to Matlab 2006b sparse matrix interface
 *  2007-02-22: Added edge centrality output
 *  2007-04-18: Updated to support additional 'weight' parameter
 *  2007-07-11: Fixed bug with weight parameter and non-zero length
 */

#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"

#include <math.h>
#include <stdlib.h>

/*
 * The mex function runs a betweenness-centrality problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int rval;
    mwIndex mrows, ncols;

    mwIndex n,nz;

    /* sparse matrix */
    mwIndex *ia, *ja;
    double *a;

    /* output data */
    double *bc;
    double *ec;

    /* used to switch between algorithm types */
    int weight_type; /* = 0 if unweighted,
                        = 1 if a vector of weights,
                        = 2 if given by the matrix */

    /*
     * The current calling pattern is
     * betweenness_centrality_mex(A,weight)
     * where weight = 0 to use the unweighted version
     *       weight = 'matrix' to use the values in the matrix
     *       weight = vector to use a vector of weights
     */

    const mxArray* arg_matrix;
    const mxArray* arg_weight;

    if (nrhs != 2)
    {
        mexErrMsgTxt("2 inputs required.");
    }

    arg_matrix = prhs[0];
    arg_weight = prhs[1];

    /* The second input must be a scalar, a string, or a vector. */
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
    if (mrows != ncols ||
        !mxIsSparse(arg_matrix) ||
        ((!mxIsDouble(arg_matrix) || mxIsComplex(arg_matrix)) && (weight_type == 2))
        )
    {
        mexErrMsgTxt("Input must be a square sparse matrix.");
    }

    n = mrows;


    /* Get the sparse matrix */

    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);

    nz = ia[n];

    /* check the reweighting array to make sure it is acceptable */
    if (weight_type == 1 && mxGetNumberOfElements(arg_weight) < nz) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexParameter",
            "the weight array must have length >= nnz(A)");
    }


    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    /* create the output vectors */
    bc = mxGetPr(plhs[0]);

    /* if they requested edge centrality, compute that as well */
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(nz,1,mxREAL);
        ec = mxGetPr(plhs[1]);
    } else {
        ec = NULL;
    }


    #ifdef _DEBUG
    mexPrintf("betweenness_centrality...");
    #endif

    rval = betweenness_centrality(n, ja, ia, a,
        bc, ec);

    #ifdef _DEBUG
    mexPrintf("done, rval=%i\n", rval);
    #endif

    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}

