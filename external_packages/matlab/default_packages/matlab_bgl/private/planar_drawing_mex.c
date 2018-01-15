/** @file planar_drawing_mex.c
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
  mwIndex nedges= 0, npos= 0;
  int rval;
  mxArray *oi, *oj, *op, *oX;
  double *i, *j, *p, *X;
  int is_maximal= 0, just_ordering= 0;

  /* [ei,ej,p,X] = planar_drawing_mex(A,is_maximal,just_ordering)
   * ei, ej is the extra set of edges
   * p is the ordering permutation
   * X is the matrix of positions
   */

  load_graph_arg(nrhs, prhs, 0, -1, -1, 0, &n, &nz, &ia, &ja, NULL);
  is_maximal= (int)load_scalar_double_arg(nrhs, prhs, 1);
  just_ordering= (int)load_scalar_double_arg(nrhs, prhs, 2);

  if (is_maximal) { nedges= 0; }
  else { nedges= n>4 ? 3*n-6 : 6; }
  plhs[0]= oi= mxCreateDoubleMatrix(nedges,1,mxREAL); i= mxGetPr(oi);
  plhs[1]= oj= mxCreateDoubleMatrix(nedges,1,mxREAL); j= mxGetPr(oj);

  plhs[2]= op= mxCreateDoubleMatrix(n,1,mxREAL); p= mxGetPr(op);
  if (just_ordering) { npos=0; } else { npos=2; }
  plhs[3]= oX= mxCreateDoubleMatrix(n,npos,mxREAL); X= mxGetPr(oX);

  nedges= 0;
  if (n==0) {
    rval= 0;
    nedges= 0;
  } else if (n==1) {
    rval= 0;
    nedges= 0;
    p[0]= 0;
    X[0]= 0; X[1]= 0;
  } else {
    rval= chrobak_payne_straight_line_drawing(n, ja, ia,
        just_ordering, is_maximal,
        (mbglIndex*)i, (mbglIndex*)j, &nedges,
        (mbglIndex*)p, (mbglDegreeType*)X);
  }

#if _DEBUG
  if (npos>0) {
    mbglDegreeType* pos= (mbglDegreeType*)X; int i, end= n>10?10:n;
    for (i=0; i<end; i++) { mexPrintf("p[%i]=(%u,%u)\n",i,pos[i],pos[i+n]); }
  }
#endif 

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
  expand_index_to_double((mwIndex*)p, p, n, 1.0);
  /* needs to be done as two calls to fix 32-bit offset errors */ 
  expand_degree_to_double((mbglDegreeType*)X, X, n*npos, 0.0);
}



