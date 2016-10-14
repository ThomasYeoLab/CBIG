/*
 * ==============================================================
 * astar_search_mex.c The mex interface to the astar_search
 * wrapper.
 *
 * David Gleich
 * 31 May 2006
 * =============================================================
 */

/*
 * 31 May 2006
 * Initial version
 *
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 *
 * 28 February 2007
 * Tried to fix bugs with Matlab 2006b sparse matrix interface.
 *
 * 19 April 2007
 * Fixed error with invalid vertex by checking the input to make sure the
 * start vertex is valid.
 *
 * Added target vertex
 *
 * 23 April 2007
 * Added edge-weight option
 */

/*
 * additional documentation for visitor for astar visitor
 *
 */

#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
#endif /* MX_API_VER */

#include "matlab_bgl.h"
#include "visitor_macros.h"
#include "expand_macros.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

PROTOTYPE_VISITOR_VERTEX_FUNCTION(initialize_vertex)
PROTOTYPE_VISITOR_VERTEX_FUNCTION(examine_vertex)
PROTOTYPE_VISITOR_VERTEX_FUNCTION(discover_vertex)
PROTOTYPE_VISITOR_VERTEX_FUNCTION(finish_vertex)

PROTOTYPE_VISITOR_EDGE_FUNCTION(examine_edge)
PROTOTYPE_VISITOR_EDGE_FUNCTION(edge_relaxed)
PROTOTYPE_VISITOR_EDGE_FUNCTION(edge_not_relaxed)
PROTOTYPE_VISITOR_EDGE_FUNCTION(black_target)


#include "expand_macros.h"

/* the heuristic function */
double astar_heuristic(void *pdata, mwIndex u);

/*
 * The mex function runs a shortest path problem.
 */
void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwIndex mrows, ncols;

  mwSize n,nz;

  /* sparse matrix */
  mwIndex *ia, *ja;
  double *a;

  /* true if this function is reweighted */
  int reweighted = 0;

  /* source */
  mbglIndex u;

  /* target */
  mwIndex v;

  double dinf;

  /* true if this function was called with a visitor */
  int use_visitor = 0;

  /* true if this function was called with a vector as h */
  int use_vector = 0;

  /* output data */
  double *d, *pred, *f;

  /*
   * The current calling pattern is
   * astar_search_mex(A,u,v,h,dinf,reweight,[visitor])
   * so visitor is an optional paramter.
   * h is either a vector with the heuristic specified for all nodes or
   *   a function handle that computes the heuristic on the fly
   * reweight is either a string of a length nnz vector
   *
   * if reweight is a length nnz vector, then we use that as the values
   * for the matrix, if its a string, then we use the values from the
   * the matrix.
   */

  const mxArray* arg_matrix;
  const mxArray* arg_source;
  const mxArray* arg_target;
  const mxArray* arg_h;
  const mxArray* arg_dinf;
  const mxArray* arg_reweight;
  const mxArray* arg_visitor = NULL;

  if (nrhs < 6 || nrhs> 7)
  {
    mexErrMsgTxt("6 or 7 inputs required.");
  }

  arg_matrix = prhs[0];
  arg_source = prhs[1];
  arg_target = prhs[2];
  arg_h = prhs[3];
  arg_dinf = prhs[4];
  arg_reweight = prhs[5];

  /* The 5th input must be a structure. */
  if (nrhs == 7)
  {
    use_visitor = 1;
    arg_visitor = prhs[6];
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

  if (use_visitor && !mxIsStruct(arg_visitor))
  {
    mexErrMsgTxt("Invalid structure.");
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

  /* Get the scalar */
  u = (mbglIndex)mxGetScalar(arg_source);
  u = u-1;

  /* Get the target */
  v = (mwIndex)mxGetScalar(arg_target);
  if (v == 0) {
    /* they want the entire search */
    v = n;
  } else {
    /* they want to look for vertex v */
    v -= 1;
  }

  if (u < 0 || u >= n)
  {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "start vertex (%i) not a valid vertex.", u+1);
  }

  if (v < 0 || v> n)
  {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "target vertex (%i) not a valid vertex.", v+1);
  }

  /* Get the uninitialized value */
  dinf = mxGetScalar(arg_dinf);

  if (mxIsDouble(arg_h))
  {
    if (mxGetNumberOfElements(arg_h) != n)
    {
      mexErrMsgTxt("The size of the vector heuristic must be the number of vertices.");
    }
    use_vector = 1;
  }
  else
  {
    int clsId = mxGetClassID(arg_h);
    if (!(clsId == mxFUNCTION_CLASS || clsId == 23))
    {
      mexErrMsgTxt("The heuristic must be either a double vector or a function.");
    }
  }

  plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,n,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);

  /* create the output vectors */

  d = mxGetPr(plhs[0]);
  pred = mxGetPr(plhs[1]);
  f = mxGetPr(plhs[2]);

#ifdef _DEBUG
  mexPrintf("astar_search...");
#endif
  if (use_visitor)
  {
    const mxArray *vis = arg_visitor;
    astar_visitor_funcs_t as_vis = {0};

    /* Check the visitor and construct the visitor structure. */
    as_vis.pdata = (void*)vis;
    CHECK_AND_SET_VISITOR_FUNCTION(vis,initialize_vertex,as_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,discover_vertex,as_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_vertex,as_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,finish_vertex,as_vis);

    CHECK_AND_SET_VISITOR_FUNCTION(vis,examine_edge,as_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_relaxed,as_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,edge_not_relaxed,as_vis);
    CHECK_AND_SET_VISITOR_FUNCTION(vis,black_target,as_vis);

    astar_search_hfunc_visitor(n, ja, ia, a,
        u,
        d, (mwIndex*)pred, f, astar_heuristic, (void*)arg_h, dinf, as_vis);
  }
  else
  {
    if (use_vector)
    {
      astar_search(n, ja, ia, a,
          u, v,
          d, (mwIndex*)pred, f, mxGetPr(arg_h), dinf);
    }
    else
    {
      astar_search_hfunc(n, ja, ia, a,
          u, v,
          d, (mwIndex*)pred, f, astar_heuristic, (void*)arg_h, dinf);
    }
  }

#ifdef _DEBUG
  mexPrintf("done\n");
#endif

  expand_index_to_double_zero_equality((mwIndex*)pred, pred, n, 1.0);

#ifdef _DEBUG
  mexPrintf("return\n");
#endif
}

double astar_heuristic(void *pdata, mwIndex u)
{
    const mxArray *f = pdata;

    mxArray* prhs[2];
    mxArray* plhs[1];

    prhs[0] = (void*)f;
    prhs[1] = mxCreateDoubleScalar((double)(u+1));

    plhs[0] = NULL;
    mexCallMATLAB(1,plhs,2,prhs,"feval");

    if (plhs[0] == NULL)
    {
        mexErrMsgTxt("Heuristic function did not return a value.");
    }

    /* mexPrintf("h(%i) = %f\n", u, (double)mxGetScalar(plhs[0])); */

    return (double)mxGetScalar(plhs[0]);
}

CALL_MATLAB_VERTEX_VISITOR_FUNCTION(initialize_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(examine_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(discover_vertex)
CALL_MATLAB_VERTEX_VISITOR_FUNCTION(finish_vertex)

CALL_MATLAB_EDGE_VISITOR_FUNCTION(examine_edge)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(edge_relaxed)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(edge_not_relaxed)
CALL_MATLAB_EDGE_VISITOR_FUNCTION(black_target)
