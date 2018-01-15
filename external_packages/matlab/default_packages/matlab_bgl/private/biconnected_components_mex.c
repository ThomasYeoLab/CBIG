/*
 * ==============================================================
 * bfs_mex.c The mex interface to the matlab bgl wrapper.
 *
 * David Gleich
 * 20 April 20020
 * =============================================================
 */

/*
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 *
 * 25 February 2007
 * Updated to use expand macros
 *
 * 1 March 2007 
 * Fixed compile bugs on win32
 */


#include "mex.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
typedef int mwSize;
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif /* MX_API_VER */

#include "matlab_bgl.h"
#include "expand_macros.h"

#include <math.h>
#include <stdlib.h>

/*
 * The mex function runs a max-flow min-cut problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int i;
    
    int mrows, ncols;
    
    int n,nz;
    
    /* sparse matrix */
    mwIndex *ia, *ja;
    
    
    /* output data */
    double *a, *ci;
    mwIndex *int_a;
     
    
    if (nrhs != 1) 
    {
        mexErrMsgTxt("1 inputs required.");
    }

    /* The first input must be a sparse matrix. */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if (mrows != ncols ||
        !mxIsSparse(prhs[0]) ||
        !mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0])) 
    {
        mexErrMsgTxt("Input must be a noncomplex square sparse matrix.");
    }
    
    n = mrows;
        
    
    
    /* Get the sparse matrix */
    
    /* recall that we've transposed the matrix */
    ja = mxGetIr(prhs[0]);
    ia = mxGetJc(prhs[0]);
    
    nz = ia[n];
        
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    a = mxGetPr(plhs[0]);
    ci = NULL;
    
    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(nz,1,mxREAL);
        ci = mxGetPr(plhs[1]);
    }
    
    int_a = (mwIndex*)a;
    for (i=0; i<n; i++)
    {
        int_a[i] = MWINDEX_MAX;
    }
    
    
    biconnected_components(n, ja, ia,
        (mwIndex*)a, (mwIndex*)ci);
    
    expand_index_to_double_zero_special((mwIndex*)a, a, n, 1.0, MWINDEX_MAX);
    if (nlhs > 1)
    {
        expand_index_to_double((mwIndex*)ci, ci, nz, 1.0);
    }
}


