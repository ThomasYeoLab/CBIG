/** @file gursoy_atun_mex.c
 * @copyright Stanford University, 2008
 * @author David F. Gleich
 * The mex interface to the matlab bgl wrapper.
 */

/** History
 * 2008-09-28: Initial version
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

/** Parse the string topology arguments
 *  [{'square'} | 'heart' | 'sphere' | 'ballN'* | 'cubeN'*]
 */
void parse_topology_arg(const char* topology,
    gursoy_atun_layout_topology_t *topo_opt, int *topo_dim)
{
  int check_strlen = 0;
  if (strstr(topology,"heart")==topology) {
    *topo_dim = 2;
    *topo_opt = HEART_LAYOUT_TOPOLOGY;
  } else if (strstr(topology,"square")==topology) {
    *topo_dim = 2;
    *topo_opt = CUBE_LAYOUT_TOPOLOGY;
  } else if (strstr(topology,"sphere")==topology) {
    *topo_dim = 3;
    *topo_opt = BALL_LAYOUT_TOPOLOGY;
  } else if (strstr(topology,"circle")==topology) {
    *topo_dim = 2;
    *topo_opt = BALL_LAYOUT_TOPOLOGY;
  } else if (strstr(topology,"ball")==topology) {
    check_strlen = 4;
    *topo_dim = atoi(topology+4);
    *topo_opt = BALL_LAYOUT_TOPOLOGY;
  } else if (strstr(topology,"cube")==topology) {
    check_strlen = 4;
    *topo_dim = atoi(topology+4);
    *topo_opt = CUBE_LAYOUT_TOPOLOGY;
  } else {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
      "the topology \'%s\' did not match the current list", topology);
  }
  if (check_strlen > 0 && strlen(topology)<=check_strlen) {
    mexErrMsgIdAndTxt("matlab_bgl:invlaidParameter",
      "the topology \'%s\' requires a dimension", topology);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwIndex n, nz, *ia, *ja; /* sparse matrix */
  double *a = NULL;
  double *X; /* output data */
  const char *topology;
  gursoy_atun_layout_topology_t topology_opt;
  double d_i=0., d_f=0., lr_i=0., lr_f=0.;
  int maxiter, progressive = 0, topology_dim;
  int rval = 0;

  /* current calling pattern:
   * gursoy_atun_mex(A,topology,maxiters,diam_i,diam_f,lr_i,lr_f,
   *   progressive_opt, edge_weights, edge_weight_opt)
   */
  const mxArray* arg_topology, *arg_maxiter;
  const mxArray* arg_diam_i, *arg_diam_f, *arg_lr_i, *arg_lr_f;
  const mxArray* arg_progressive_opt;
  if (nrhs != 10) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall","10 inputs required.");
  }
  arg_topology= prhs[1]; arg_maxiter= prhs[2];
  arg_diam_i= prhs[3]; arg_diam_f= prhs[4];
  arg_lr_i= prhs[5]; arg_lr_f= prhs[6];
  arg_progressive_opt= prhs[7];

  load_graph_arg(nrhs, prhs, 0, 8, 9, 1, &n, &nz, &ia, &ja, &a);

  /* get the parameters */
  topology= load_string_arg(arg_topology,1);
  parse_topology_arg(topology, &topology_opt, &topology_dim);

  if (isscalardouble(arg_maxiter)
      && isscalardouble(arg_diam_i) && isscalardouble(arg_diam_f)
      && isscalardouble(arg_lr_i) && isscalardouble(arg_lr_f)) {
    maxiter = mxGetScalar(arg_maxiter);
    d_i = mxGetScalar(arg_diam_i); d_f = mxGetScalar(arg_diam_f);
    lr_i = mxGetScalar(arg_lr_i); lr_f = mxGetScalar(arg_lr_f);
  } else {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The scalar parameters must be scalars of type double");
    return;
  }
  /* check if they want to reuse old positions */
  if (!mxIsEmpty(arg_progressive_opt)) {
    if (mxGetM(arg_progressive_opt) != n
        || mxGetN(arg_progressive_opt) != topology_dim) {
      mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
          "The progressive call requires an %i-by-%i input matrix of positions",
          n, topology_dim);
    }
    plhs[0]= mxDuplicateArray(arg_progressive_opt);
    progressive = 1;
  } else {
    plhs[0]= mxCreateDoubleMatrix(n,topology_dim,mxREAL);
  }
  X= mxGetPr(plhs[0]);
  if (n==1) { return; } /* special case singleton graph */
  rval= gursoy_atun_layout(n, ja, ia, a, topology_opt, topology_dim,
    maxiter, d_i, d_f, lr_i, lr_f, X);
  if (rval == gursoy_atun_invalid_dim) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
      "The topology dimension %i was larger than the precompiled code allows",
      topology_dim);
  } else if (rval == gursoy_atun_invalid_dim) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
      "The topology dimension %i is not valid (maybe too small?)",
      topology_dim);
  } else if (rval != 0) {
    mexErrMsgIdAndTxt("matlab_bgl:callFailed",
        "The libmbgl call failed with rval=%i", rval);
  }
}
