#ifndef MATLAB_BGL_COMMON_FUNCTIONS_H
#define MATLAB_BGL_COMMON_FUNCTIONS_H

/** @file common_functions.h
 * @copyright Stanford University, 2008
 * @author David F. Gleich
 * Implement a few common functions to be included in many mex files.
 */

/** History
 *  2007-07-08: Initial version
 *  2008-09-26: Added isscalar and isscalardoube
 *  2008-09-27: Added intmin
 */

int intmin(int a, int b) { if (a<b) { return a; } else { return b; } }

int isscalar(const mxArray* a) {
	if (mxGetM(a) == 1 && mxGetN(a) == 1 && mxGetNumberOfDimensions(a) == 2) {
		return 1;
	} else {
		return 0;
	}
}

int isscalardouble(const mxArray* a) {
	if (isscalar(a) && mxIsDouble(a)) {
		return 1;
	} else {
		return 0;
	}
}

char* load_string_arg(const mxArray* a, int k)
{
    mwSize buflen;
    char *s;
    int status;

    if (mxIsChar(a) != 1) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "argument %i must be a string", k+1);
    }

    /* Input must be a row vector. */
    if (mxGetM(a) != 1) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "argument %i must be a string (and a row vector)", k+1);
    }

    /* Get the length of the input string. */
    buflen = (mxGetM(a) * mxGetN(a)) + 1;

    /* Allocate memory for input and output strings. */
    s = mxCalloc(buflen, sizeof(char));

    status = mxGetString(a, s, buflen);
    if (status != 0) {
        mexErrMsgIdAndTxt("matlab_bgl:sizeError",
            "insufficient space to copy argument %i to a string", k+1);
    }

    return s;
}

double load_scalar_arg(const mxArray* a, int k)
{
    if (mxGetNumberOfElements(a) > 1 || !mxIsDouble(a)) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "argument %i should be a scalar", k+1);
    }
    return mxGetScalar(a);
}

double load_scalar_double_arg(int nrhs, const mxArray *prhs[], int arg) {
  const mxArray *a;
  if (nrhs <= arg) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall",
        "load_scalar on arg %i failed because only %i args given",
        arg, nrhs);
  }
  a= prhs[arg];
  if (mxGetNumberOfElements(a) != 1 || !mxIsDouble(a)) {
      mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
          "argument %i should be a scalar", arg+1);
  }
  return mxGetScalar(a);
}

double* load_vector_double_arg(int nrhs, const mxArray *prhs[], int arg,
    mwIndex minlen, mwIndex maxlen, int check_min, int check_max) {
  const mxArray *a;
  mwSize numel;
  if (nrhs <= arg) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall",
        "load_vector on arg %i failed because only %i args given",
        arg, nrhs);
  }
  a= prhs[arg];
  if (!mxIsDouble(a) || mxIsComplex(a)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "vector arg=%i must have type double.", arg);
  }
  numel = mxGetNumberOfElements(a);
  if (check_min && numel<minlen) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
      "vector arg=%i must have numel() at least %i", arg, minlen);
  }
  if (check_max && numel>maxlen) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
      "vector arg=%i must have numel() at most %i", arg, maxlen);
  }
  return mxGetPr(a);
}

double* load_matrix_double_arg(int nrhs, const mxArray *prhs[],
    int arg, const char* argname,
    int empty_okay, mwSize *mrows, mwSize *ncols,
    int exact_size, mwSize m, mwSize n)
{
  const mxArray *a;
  if (nrhs <= arg) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall",
        "load_matrix on %s (arg=%i) failed because only %i args given",
        argname, arg, nrhs);
  }
  a= prhs[arg];
  if (mxGetNumberOfDimensions(a) != 2) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
      "matrix %s (arg=%i) must have only 2 dimension.", argname, arg);
  }
  if (!mxIsDouble(a) || mxIsComplex(a)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
      "matrix %s (arg=%i) must have type double.", argname, arg);
  }
  if (exact_size && m!=mxGetM(a) && n!=mxGetN(a)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
      "matrix %s (arg=%i) must have size %i-by-%i.", argname, arg, m, n);
  }
  if (!empty_okay && mxIsEmpty(a)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
     "matrix %s (arg=%i) cannot be empty.", argname, arg);
  }
  if (mrows) { *mrows = mxGetM(a); }
  if (ncols) { *ncols = mxGetN(a); }
  return mxGetPr(a);
}

void load_graph_arg(int nrhs, const mxArray *prhs[],
  int arg, int reweighted_arg, int reweighted_opt_arg,
  int needs_weights,
  mwIndex *nverts, mwIndex *nedges,
  mwIndex **ia, mwIndex **ja, double **a)
{
  mwIndex m, n, *ir, *jc;
  const mxArray* arg_matrix;
  int reweighted = 0;

  if (nrhs <= arg || nrhs <= reweighted_arg || nrhs <= reweighted_opt_arg) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidCall",
        "load_matrix requires args %i,%i,%i but only %i args provided",
        arg, reweighted_opt_arg, reweighted_arg, nrhs);
  }

  arg_matrix= prhs[arg];
  if (needs_weights && reweighted_arg>=0 && reweighted_opt_arg>=0) {
    reweighted = (int)load_scalar_double_arg(nrhs,prhs,reweighted_arg);
  }
  m = mxGetM(arg_matrix);
  n = mxGetN(arg_matrix);
  if (m != n || !mxIsSparse(arg_matrix)) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The graph (arg=%i) must be square and sparse.", arg);
  }
  if (needs_weights &&
      (!reweighted && (!mxIsDouble(arg_matrix) || mxIsComplex(arg_matrix)))) {
    mexErrMsgIdAndTxt("matlab_bgl:invalidParameter",
        "The graph (arg=%i) must have type double.", arg);
  }

  ir = mxGetIr(arg_matrix);
  jc = mxGetJc(arg_matrix);

  if (nverts) { *nverts = n; }
  if (nedges) { *nedges = jc[n]; }
  if (ja) { *ja = ir; }
  if (ia) { *ia = jc; }

  if (needs_weights && a) {
    if (reweighted) {
      *a = load_vector_double_arg(nrhs, prhs, reweighted_opt_arg, jc[n], 0, 1, 0);
    } else {
      *a = mxGetPr(arg_matrix);
    }
  } else {
    if (a) { *a = NULL; }
  }
}


#endif /* MATLAB_BGL_COMMON_FUNCTIONS_H */

