/*
 * ==============================================================
 * path_from_pred.c A mex interface to the matlab bgl function
 * to convert a predecessor array into a path
 *
 * David Gleich
 * 17 April 2007
 * =============================================================
 */

/*
 * History
 *
 * 17 April 2007
 * Initial version
 */


#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include <math.h>

/*
 * The mex function traces back a path one vertex to a predecessor vertex.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize mrows, ncols;

    /* size of pred array */
    mwSize n;
    double *pred;

    /* dest */
    mwIndex d;

    /* output data */
    mwSize pathlen;
    double *path;

    if (nrhs != 2)
    {
        mexErrMsgTxt("2 inputs required.");
    }

    /* The first input must be an array. */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if ((mrows != 1 && ncols != 1) ||
        !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("Predecessor input must be a double array.");
    }

    n = mxGetNumberOfElements(prhs[0]);

    /* The second input must be a scalar. */
    if (mxGetNumberOfElements(prhs[1]) > 1 || !mxIsDouble(prhs[1]))
    {
        mexErrMsgTxt("Invalid scalar for destination.");
    }

    pred = mxGetPr(prhs[0]);
    d = (mwIndex)mxGetScalar(prhs[1]);

    /* Adjust the matlab index to a c index */
    d = d - 1;

    /* Compute the path length and verify the predecessory array */
    pathlen = 0;
    while (d >= 0 && d < n && pathlen < n) {
         d = (mwIndex)floor(pred[d]);
         pathlen++;

         if (d == 0) {
             /* found the source */
             break;
         }
         d -= 1;
    }

    if (d != 0) {
        /* this case indicates that the predecessor array was incorrectly
         * formed, so throw an error */
        if (d >= n) { mexErrMsgTxt("The predecessor array referenced a vertex too large."); }
        else if (d < -1) { mexErrMsgTxt("The predecessor array contains negative references."); }
        else if (pathlen > n-1) { mexErrMsgTxt("The predecessor array contained a loop."); }
        else { mexErrMsgTxt("Unknown error with predecessor array."); }
    }

    /* Now allocate and build the path */
    plhs[0] = mxCreateDoubleMatrix(1,pathlen,mxREAL);
    path = mxGetPr(plhs[0]);

    /* Get the ending index again */
    d = (mwIndex)mxGetScalar(prhs[1]);

    /* Adjust the matlab index to a c index */
    d = d - 1;

    while (d >= 0 && d < n && pathlen > 0 && pathlen <= n) {
        path[pathlen-1] = (double)(d+1);
        d = (mwIndex)floor(pred[d]);
        pathlen--;

        if (d == 0) {
            /* found the source */
            break;
        }

        d -= 1;
    }

    if (pathlen != 0) { mexErrMsgTxt("Unknown error creating path."); }
}

