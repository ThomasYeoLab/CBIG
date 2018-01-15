/** @file kamada_kawai_spring_layout_mex.c
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
  mwIndex mrows, ncols, n, nz;
  mwIndex *ia, *ja; /* sparse matrix */
  double *a = NULL;
  int reweighted = 0; /* true if this function is reweighted */
  double *X, *S, *D; /* output data */
  double tol, k, edgelen;
  int maxiter= 100, progressive = 0;
  int rval = 0;

  /* current calling pattern:
   * kamada_kawai_spring_layout_mex(A,tol,maxiter,spring_constant,...
   *   progressive_opt,edge_length, edge_weights, edge_weight_opt)
   */
  const mxArray* arg_matrix;
  const mxArray* arg_tol, *arg_maxiter, *arg_k, 
          *arg_progressive_opt, *arg_edgelen;
  const mxArray* arg_ews, *arg_ewopt;
  if (nrhs != 8) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall","8 inputs required.");
  }
  arg_matrix= prhs[0];
  arg_tol= prhs[1]; arg_maxiter= prhs[2]; arg_k= prhs[3];
  arg_progressive_opt= prhs[4]; arg_edgelen= prhs[5];
  arg_ews= prhs[6]; arg_ewopt=prhs[7];
  if ((int)mxGetScalar(arg_ews)!=0) { reweighted=1; }
  /* The first input must be a sparse matrix. */
  mrows = mxGetM(arg_matrix);
  ncols = mxGetN(arg_matrix);
  if (mrows != ncols || !mxIsSparse(arg_matrix)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The graph must be square and sparse.");
  }
  if (!reweighted && (!mxIsDouble(arg_matrix) || mxIsComplex(arg_matrix))) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The graph must have type double.");
  }

  /* Get the sparse matrix, and recall that we've transposed the matrix */
  n = mrows;
  ja = mxGetIr(arg_matrix);
  ia = mxGetJc(arg_matrix);
  nz = ia[n];
  if ((int)mxGetScalar(arg_ews)==0) { a = mxGetPr(arg_matrix); }
  else if ((int)mxGetScalar(arg_ews)==1) {
    if (mxGetNumberOfElements(arg_ewopt) < nz) {
      mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
          "The reweight array must have length at least nnz(A)");
    }
    if (!mxIsDouble(arg_ewopt) < nz) {
      mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
          "The reweight array must be of type double");
    }
    a= mxGetPr(arg_ewopt);
  }
  /* get the parameters */
  if (isscalardouble(arg_tol) && isscalardouble(arg_maxiter) 
      && isscalardouble(arg_k) && isscalardouble(arg_edgelen)) {
    tol = mxGetScalar(arg_tol);
    maxiter = mxGetScalar(arg_maxiter);
    k = mxGetScalar(arg_k);
    edgelen = mxGetScalar(arg_edgelen);
  } else {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The scalar parameters must be scalars of type double");
    return;
  }
  /* create the output vectors */
  plhs[1]= mxCreateDoubleMatrix(n,n,mxREAL);
  plhs[2]= mxCreateDoubleMatrix(n,n,mxREAL);
  S= mxGetPr(plhs[1]);
  D= mxGetPr(plhs[2]);
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
  if (n==0) { return; } /* special case empty graph */
  if (n==1) { return; } /* special case singleton graph */
  rval= kamada_kawai_spring_layout(n, ja, ia, a, tol, maxiter, k, progressive,
           edgelen, X, S, D);
  if (rval == -2) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
      "the graph must be connected and have no negative weight cycles");
  } else if (rval != 0) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
        "The libmbgl call failed with rval=%i", rval);
  }
}
