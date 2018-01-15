/** @file planar_test_mex.c
 * @copyright Stanford University, 2008
 * @author David F. Gleich
 * The mex interface to the matlab bgl wrapper.
 */

/** History
 * 2008-10-01: Initial version
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
  mwIndex n, nz, *ia, *ja; /* sparse matrix */
  int test_type = 0, rval=-1;

  load_graph_arg(nrhs, prhs, 0, -1, -1, 0, &n, &nz, &ia, &ja, NULL);
  test_type = (int)load_scalar_double_arg(nrhs, prhs, 1);

  /* [is_planar ksubgraph embedding] = planar_test_mex(A,0)
   * [is_kuratowski] = planar_test_mex(A,1)
   * [is_straight_line] = planar_test_mex(A,2,X)
   */

  plhs[0]= mxCreateDoubleMatrix(1,1,mxREAL);

  if (test_type == 0)
  {
    /* get planar subgraph information */
    int is_planar = 0;

    if (nlhs <= 1) {
      /* just test for the planar graph */
      rval= boyer_myrvold_planarity_test(n, ja, ia, &is_planar,
          NULL, NULL, NULL, NULL, NULL);
    } else if (nlhs <= 3) {
      /* test for the planar graph and the kuratowski subgraph */
      double *ki, *kj;
      mwIndex nedges= n>4 ? 3*n-6 : 6;
      plhs[1]= mxCreateDoubleMatrix(nedges,1,mxREAL); ki=mxGetPr(plhs[1]);
      plhs[2]= mxCreateDoubleMatrix(nedges,1,mxREAL); kj=mxGetPr(plhs[2]);
      nedges= 0;
      rval= boyer_myrvold_planarity_test(n, ja, ia, &is_planar,
          (mwIndex*)ki, (mwIndex*)kj, &nedges, NULL, NULL);
      expand_index_to_double((mwIndex*)ki, ki, nedges, 1.0);
      expand_index_to_double((mwIndex*)kj, kj, nedges, 1.0);
      mxSetM(plhs[1], nedges);
      mxSetM(plhs[2], nedges);
    } else if (nlhs > 3) {
      /* test for the planar graph and the edge order */
      double *eip, *eie;
      double *ki, *kj;
      mwIndex nedges= n>4 ? 3*n-6 : 6;
      plhs[1]= mxCreateDoubleMatrix(nedges,1,mxREAL); ki=mxGetPr(plhs[1]);
      plhs[2]= mxCreateDoubleMatrix(nedges,1,mxREAL); kj=mxGetPr(plhs[2]);
      nedges= 0;
      plhs[3]= mxCreateDoubleMatrix(n+1,1,mxREAL); eip= mxGetPr(plhs[3]);
      plhs[4]= mxCreateDoubleMatrix(nz,1,mxREAL); eie= mxGetPr(plhs[4]);
      rval= boyer_myrvold_planarity_test(n, ja, ia, &is_planar,
          (mwIndex*)ki, (mwIndex*)kj, &nedges, (mwIndex*)eip, (mwIndex*)eie);
      expand_index_to_double((mwIndex*)ki, ki, nedges, 1.0);
      expand_index_to_double((mwIndex*)kj, kj, nedges, 1.0);
      mxSetM(plhs[1], nedges);
      mxSetM(plhs[2], nedges);
      expand_index_to_double((mwIndex*)eip, eip, n+1, 1.0);
      expand_index_to_double((mwIndex*)eie, eie, nz, 1.0);
    }
    expand_int_to_double(&is_planar, mxGetPr(plhs[0]), 1, 0.0);
  }
  else if (test_type == 1)
  {
    /* test for a kuratowski subgraph */
    int is_ksubgraph= 0;
    rval= is_kuratowski_subgraph(n, ja, ia, &is_ksubgraph);
    expand_int_to_double(&is_ksubgraph, mxGetPr(plhs[0]), 1, 0.0);
  }
  else if (test_type == 2)
  {
    /** Test for a straight line embedding */
    int is_sldrawing= 0;
    double* X= load_matrix_double_arg(nrhs,prhs,2,"X",1,NULL,NULL,1,n,2);
    rval= is_straight_line_drawing(n, ja, ia, X, &is_sldrawing);
    if (rval==-11) {
      mexErrMsgIdAndTxt("matlab_bgl:callFailed",
          "is_straight_line_drawing requires positive positions.");
    }
    expand_int_to_double(&is_sldrawing, mxGetPr(plhs[0]), 1, 0.0);
  }

  if (rval != 0) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
        "The libmbgl call failed with rval=%i", rval);
  }
}
