/*
 * David Gleich
 * Copyright, Stanford Unviersity, 2007
 */

/** 
 * @file test_matching_mex.c
 * Wrap a call to the libmbgl test_matching function.
 */

/*
 * 8 July 2007
 * Initial version
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

/*
 * The mex function runs a MST problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{ 
    mwIndex i;
    
    mwIndex mrows, ncols;
    
    mwIndex n,nz;
    
    /* sparse matrix */
    mwIndex *ia, *ja;
    
    /* matching */
    mwIndex *m;
    double *m_double;
    
    /* output */
    int verify = 0;
    
    /* 
     * The current calling pattern is
     * test_matching_mex(A,matching)
     */
    
    const mxArray* arg_matrix;
    const mxArray* arg_matching;
    int required_arguments = 2;
    
    if (nrhs != required_arguments) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the function requires %i arguments, not %i\n", 
            required_arguments, nrhs);
    }
    
    arg_matrix = prhs[0];
    arg_matching = prhs[1];
    
    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    if (mrows != ncols || !mxIsSparse(arg_matrix)) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the matrix must be sparse and square");
    }

    n = mrows;
    
    /* The second input must be of size n */
    if (mxGetNumberOfElements(arg_matching) != n) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the matching must be size %i not %i", n,
            mxGetNumberOfElements(arg_matching));
    }
    m_double = mxGetPr(arg_matching);
    m = mxCalloc(n, sizeof(mwIndex));
    for (i=0; i < n; i++) {
        m[i] = (mwIndex)m_double[i] ;
        if (m[i] == 0) { m[i] = n; }
        else { --m[i]; }
    }
    
    /* Get the sparse matrix */
    
    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);
    
    nz = ia[n];
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

    #ifdef _DEBUG
    mexPrintf("test_maximum_cardinality_matching...");
    #endif 

    test_maximum_cardinality_matching(n, ja, ia, m, &verify);
    
    *mxGetPr(plhs[0]) = (double)verify;
    
    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif 
    
}

