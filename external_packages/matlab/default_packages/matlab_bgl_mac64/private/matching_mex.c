/*
 * David Gleich
 * Copyright, Stanford Unviersity, 2007
 */

/** 
 * @file matching_mex.c
 * Wrap a call to the libmbgl matching function.
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
    mwIndex mrows, ncols;
    
    mwIndex n,nz;
    
    /* sparse matrix */
    mwIndex *ia, *ja;
    
    int verify;
    char *initial_match_name;
    char *augmenting_path_alg_name;
    
    int initial_match = 0;
    int augmenting_path = 0;
    
    /* 
     * The current calling pattern is
     * matching_mex(A,verify,initial_match_name,augmenting_path_name)
     */
    
    const mxArray* arg_matrix;
    const mxArray* arg_verify;
    const mxArray* arg_initial_match_name;
    const mxArray* arg_augmenting_path_name;    
    int required_arguments = 4;
    
    if (nrhs != required_arguments) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the function requires %i arguments, not %i\n", 
            required_arguments, nrhs);
    }
    
    arg_matrix = prhs[0];
    arg_verify = prhs[1];
    arg_initial_match_name = prhs[2];
    arg_augmenting_path_name = prhs[3];
    
    verify = (int)mxGetScalar(arg_verify) + 1;
    initial_match_name = load_string_arg(arg_initial_match_name,2);
    augmenting_path_alg_name = load_string_arg(arg_augmenting_path_name,3);

    /* The first input must be a sparse matrix. */
    mrows = mxGetM(arg_matrix);
    ncols = mxGetN(arg_matrix);
    if (mrows != ncols || !mxIsSparse(arg_matrix)) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the matrix must be sparse and square");
    }

    n = mrows;
         
    
    /* Get the sparse matrix */
    
    /* recall that we've transposed the matrix */
    ja = mxGetIr(arg_matrix);
    ia = mxGetJc(arg_matrix);
    
    nz = ia[n];
    
    
    
  	/* convert options into numeric values */
    if (strcmp(initial_match_name,"none") == 0) {
        initial_match = 1;
    } else if (strcmp(initial_match_name,"greedy") == 0) {
        initial_match = 2;
    } else if (strcmp(initial_match_name,"extra_greedy") == 0) {
        initial_match = 3;
    } else {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "initial_match option %s is invalid\n", initial_match_name);
    }

    if (strcmp(augmenting_path_alg_name,"none") == 0) {
        augmenting_path = 1;
    } else if (strcmp(augmenting_path_alg_name,"edmonds") == 0) {
        augmenting_path = 2;
    } else {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "augmenting_path_alg option %s is invalid\n", 
            augmenting_path_alg_name);
    }
    
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

    {
        mwIndex null_vertex;
        int verified=0;
        
        #ifdef _DEBUG
        mexPrintf("matching(%s,%s,%i)...", 
            initial_match_name, augmenting_path_alg_name,verify);
        #endif 
    
        maximum_cardinality_matching(n, ja, ia, 
            (mwIndex*)mxGetPr(plhs[0]),
            initial_match, augmenting_path, verify,
            &verified, &null_vertex);
        
        #ifdef _DEBUG
        mexPrintf("done!\n");
        #endif 
     
        expand_index_to_double_zero_special(
            (mwIndex*)mxGetPr(plhs[0]), mxGetPr(plhs[0]),
            n, 1.0, null_vertex);
        
        *mxGetPr(plhs[1]) = (double)verified;
    }
    
    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif 
    
}

