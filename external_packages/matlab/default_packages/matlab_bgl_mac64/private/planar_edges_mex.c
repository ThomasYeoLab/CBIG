/** @file planar_edges_mex.c
 * @copyright Stanford University, 2008
 * @author David F. Gleich
 * Get extra sets of edges for a planar graph
 */

/** History
 * 2008-10-05: Initial version
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
  mwIndex nedges=0;
  int rval;
  mxArray *oi, *oj;
  double *i, *j;
  int make_connected, make_maximal, make_biconnected;

  /* [i j]=planar_edges_mex(G,make_connected,make_biconnected,make_maximal)
   */
  load_graph_arg(nrhs, prhs, 0, -1, -1, 0, &n, &nz, &ia, &ja, NULL);
  make_connected = (int)load_scalar_double_arg(nrhs, prhs, 1);
  make_biconnected = (int)load_scalar_double_arg(nrhs, prhs, 2);
  make_maximal = (int)load_scalar_double_arg(nrhs, prhs, 3);

  nedges= n>4 ? 3*n-6 : 6;
  plhs[0]= oi= mxCreateDoubleMatrix(nedges,1,mxREAL); i= mxGetPr(oi);
  plhs[1]= oj= mxCreateDoubleMatrix(nedges,1,mxREAL); j= mxGetPr(oj);

  nedges= 0;
  rval= triangulate_graph(n, ja, ia, make_connected, make_biconnected, make_maximal,
      (mbglIndex*)i, (mbglIndex*)j, &nedges );
  if (rval == 1) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
      "triangulate_graph with biconnected or maximal requires a planar graph");
  } else if (rval != 0) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
        "The libmbgl call failed with rval=%i", rval);
  }
  expand_index_to_double((mwIndex*)i, i, nedges, 1.0);
  expand_index_to_double((mwIndex*)j, j, nedges, 1.0);
  mxSetM(oi, nedges);
  mxSetM(oj, nedges);
}
