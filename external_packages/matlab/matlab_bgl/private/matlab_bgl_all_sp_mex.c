/*
 * ==============================================================
 * matlab_bgl_all_sp_mex.c The mex interface to the matlab bgl wrapper.
 *
 * David Gleich
 * 20 April 2006
 * =============================================================
 */

/*
 * 20 April 2006
 * Initial version
 *
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 *
 * 1 March 2007
 * Updated to get predecessors from floyd warshall and use expand macros
 *
 * 18 April 2007
 * Updated to support additional 'weight' parameter
 *
 * 2008-04-02: Fixed bug with predecessor return and off-by-1 on diagonal
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
    double *a;

    /* true if this function is reweighted */
    int reweighted = 0;

    double dinf;

    /* output data */
    double *D;
    int rval;

    /* sp type string */
    int buflen;
    char *algname;
    int status;

    /*
     * The current calling pattern is
     * matlab_bgl_all_sp_mex(A,algname,dinf,reweight)
     * algname is a string with either 'johnson', 'floyd_warshall'
     * reweight is either a string or a length nnz vector
     *
     * if reweight is a length nnz vector, then we use that as the values
     * for the matrix, if its a string, then we use the values from the
     * the matrix.
     */

    const mxArray* arg_matrix;
    const mxArray* arg_algname;
    const mxArray* arg_dinf;
    const mxArray* arg_reweight;

    if (nrhs != 4)
    {
        mexErrMsgTxt("4 inputs required.");
    }

    arg_matrix = prhs[0];
    arg_algname = prhs[1];
    arg_dinf = prhs[2];
    arg_reweight = prhs[3];

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

    /* Get the uninitialized value */
    dinf = mxGetScalar(arg_dinf);

    /* Get the algorithm type */

    if (mxIsChar(arg_algname) != 1)
        mexErrMsgTxt("Input 2 must be a string (algname).");

    /* Input must be a row vector. */
    if (mxGetM(arg_algname) != 1)
        mexErrMsgTxt("Input 2 must be a row vector.");

    /* Get the length of the input string. */
    buflen = (mxGetM(arg_algname) * mxGetN(arg_algname)) + 1;

    /* Allocate memory for input and output strings. */
    algname = mxCalloc(buflen, sizeof(char));

    status = mxGetString(arg_algname, algname, buflen);
    if (status != 0)
        mexErrMsgTxt("Not enough space for algname input.");

    plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);

    /* create the output vectors */

    D = mxGetPr(plhs[0]);


    #ifdef _DEBUG
    mexPrintf("all_sp_%s...",algname);
    #endif
    if (strcmp(algname, "johnson") == 0)
    {
        rval = johnson_all_sp(n, ja, ia, a,
            D, dinf);
    }
    else if (strcmp(algname, "floyd_warshall") == 0)
    {
        double *pred = NULL;
        if (nlhs > 1) {
            plhs[1] = mxCreateDoubleMatrix(n,n,mxREAL);
            pred = mxGetPr(plhs[1]);
        }
        rval = floyd_warshall_all_sp(n, ja, ia, a,
            D, dinf, (mbglIndex*)pred);
        if (pred) {
            mwIndex i;
            expand_index_to_double((mwIndex*)pred, pred, n*n, 1.0);
            /* zero out entries on the diagonal */
            for (i=0; i<n; i++) {
                pred[i+i*n]=0.0;
            }
        }
    }
    else
    {
        mexErrMsgTxt("Unknown algname.");
        return;
    }
    #ifdef _DEBUG
    mexPrintf("done!\n");
    #endif

    if (rval != 0)
    {
        mexErrMsgTxt("Negative weight cycle detected, check the input.");
    }

    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif
}

