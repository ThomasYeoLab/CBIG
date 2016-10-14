/** @file fruchterman_reingold_mex.c
 * @copyright Stanford University, 2008
 * @author David F. Gleich
 * The mex interface to the matlab bgl wrapper.
 */

/** History
 * 2008-09-26: Initial version
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwIndex mrows, ncols, n;
  mwIndex *ia, *ja; /* sparse matrix */
  double *X;
  int maxiter, fptype;
  double itemp, width, height;
  int progressive = 0;
  int rval;

  /* current calling pattern:
   * fruchterman_reingold_mex(A, maxiter, fptype, width, height, progressive_opt)
   */
  const mxArray* arg_matrix;
  const mxArray* arg_maxiter,*arg_itemp, *arg_fptype, *arg_width, *arg_height;
  const mxArray* arg_progressive_opt;
  if (nrhs != 7) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall","7 inputs required.");
  }
  arg_matrix= prhs[0];
  arg_maxiter= prhs[1]; arg_itemp= prhs[2]; arg_fptype= prhs[3];
  arg_width= prhs[4]; arg_height= prhs[5];
  arg_progressive_opt= prhs[6];

  /* The first input must be a sparse matrix. */
  mrows = mxGetM(arg_matrix);
  ncols = mxGetN(arg_matrix);
  if (mrows != ncols || !mxIsSparse(arg_matrix)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The graph must be square and sparse.");
  }

  /* Get the sparse matrix, and recall that we've transposed the matrix */
  n = mrows;
  ja = mxGetIr(arg_matrix);
  ia = mxGetJc(arg_matrix);
  /* get the parameters */
  if (isscalardouble(arg_maxiter) && isscalardouble(arg_itemp)
      && isscalardouble(arg_fptype) && isscalardouble(arg_width)
      && isscalardouble(arg_height)) {
    maxiter = (int)mxGetScalar(arg_maxiter);
    itemp = mxGetScalar(arg_itemp);
    fptype = (int)mxGetScalar(arg_fptype);
    width = mxGetScalar(arg_width);
    height = mxGetScalar(arg_height);
  } else {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The scalar parameters must be scalars of type double");
    return;
  }

  /* check if they want to reuse old positions */
  if (!mxIsEmpty(arg_progressive_opt)) {
    if (mxGetM(arg_progressive_opt) != n || mxGetN(arg_progressive_opt) != 2) {
      mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
          "The progressive call requires an %i-by-2 input matrix of positions",
          n);
    }
    plhs[0]= mxDuplicateArray(arg_progressive_opt);
    progressive = 1;
  } else {
    plhs[0]= mxCreateDoubleMatrix(n,2,mxREAL);
  }
  X= mxGetPr(plhs[0]);
#ifdef _DEBUG
  mexPrintf("fruchterman_reingold_force_directed_layout(%i,%i,%lf,%lf,%i)\n",
      maxiter, fptype, width, height, progressive);
  for (rval=0; rval<intmin(5,n); rval++) {
    mexPrintf("p[%i]=(%lf,%lf)\n",rval,X[rval],X[rval+n]); }
#endif /* _DEBUG */
  if (fptype == 0) {
    rval = fruchterman_reingold_force_directed_layout(n, ja, ia,
        maxiter, itemp, 0, width, height, progressive, X);
  } else {
    rval = fruchterman_reingold_force_directed_layout(n, ja, ia,
        maxiter, itemp, 1, width, height, progressive, X);
  }
#ifdef _DEBUG
  mexPrintf("return, rval = %i\n", rval);
  for (rval=0; rval<intmin(5,n); rval++) {
      mexPrintf("p[%i]=(%lf,%lf)\n",rval,X[rval],X[rval+n]); }
#endif /* _DEBUG */
}

