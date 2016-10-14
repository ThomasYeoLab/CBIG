/*
 * ==============================================================
 * max_flow_mex.c The mex interface to the matlab bgl wrapper.
 *
 * David Gleich
 * 16 April 2006
 * =============================================================
 */

/*
 * 15 June 2006 
 * Changed error message.
 *
 * 19 February 2007
 * Updated to use Matlab 2006b sparse matrix interface
 *
 * 8 July 2007
 * Added algname option and kolmogorov and edmunds interfaces
 *
 * 2008-04-02: Fix 1 for max-flow cut bug
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


/**
 * The output matrix B has twice the number of non-zeros as A, it also 
 * row oriented, not column oriented.
 *
 * For each non-zero (i,j), we insert (i,j) with capacity a(i,j) and
 * (j,i) with capacity 0.
 *
 * @param n number of nodes in matrix a
 * @param pja the column pointer (length n+1)
 * @param ia the row array (length pja[n])
 * @param a the value array (length pja[n]), we floor these to integers in b
 * @param pib the row pointer for the new matrix b (length n+1)
 * @param jb the column array (length pib[n])
 * @param cb the capacity array (length pib[n])
 * @param rev_map the reverse edge map (length pib[n])
 */
void build_matrix(mbglIndex n, mbglIndex *ia, mbglIndex *pja, double *a, mbglIndex **opib, mbglIndex **ojb, int **ocb,
            mbglIndex **orev_map)
{
    mbglIndex i,j,k;
    mbglIndex *pib, *jb, *rev_map;
    int *cb;
    
    mbglIndex nz = pja[n];
    
    /* allocate the row pointer for the new matrix */
    pib = mxCalloc(sizeof(mbglIndex),(n+1));
    jb = mxCalloc(sizeof(mbglIndex),(2*nz));
    cb = mxCalloc(sizeof(int),(2*nz));
    rev_map = mxCalloc(sizeof(mbglIndex),(2*nz));   
    
    *opib = pib;
    *ojb = jb;
    *ocb = cb;
    *orev_map = rev_map;
    
    for (j=0;j<n;j++)
    {
        for (k=pja[j];k<pja[j+1];k++)
        {
            i = ia[k];
            
            /* add the degree for each entry, */
            pib[i+1]++;
            pib[j+1]++;
        }
    }
    
    
    
    /* compute the starting locations for each row */
    k=0;
    for (i=1;i<n+1;i++)
    {
        k += pib[i];
        pib[i] = k;
    }
    
    
    /* quick check... */
    if (pib[n] != 2*nz)
    {
        mexErrMsgTxt("Error: number of non-zeros do not match");
    }
    
    

    /* now actually add all the data to the matrix */
    for (j=0;j<n;j++)
    {
        for (k=pja[j];k<pja[j+1];k++)
        {
            i = ia[k];
            
            /* insert edge (i,j) */
            jb[pib[i]] = j;
            cb[pib[i]] = (int)floor(a[k]);
            
            /* insert edge (j,i) */
            jb[pib[j]] = i;
            cb[pib[j]] = 0;
            
            /* create reverse map */
            rev_map[pib[i]] = pib[j];
            rev_map[pib[j]] = pib[i];
            
            /* increment pointers */
            pib[i]++;
            pib[j]++;
        }
    }
    
    
    /* now restore the correct value on the array, that is, shift it by one */
    for (i=n;i>0;i--)
    {
        pib[i] = pib[i-1];
    }
    
    pib[0] = 0;
}

/**
 * Convert the residual graph into a graph cut.
 * 
 * @param u the source vertex
 * @param n the number of nodes
 * @param pia the row pointer for the flow graph A
 * @param ja the column array for flow graph A
 * @param cap the capacity array for flow graph A
 * @param res the residual array for flow graph A
 * @param mincut the mincut array (output)
 */
void build_cut(mbglIndex u, mbglIndex n, mbglIndex *pia, mbglIndex *ja, int *cap, int *res, int *mincut)
{
    mbglIndex k;
    mbglIndex v, vn;
    
    mbglIndex *q;
    mbglIndex qhead, qtail;
    
    /* initialize the mincut array to 0 */
    memset(mincut, 0, sizeof(int)*n);
    
    /* 
     * if mincut[i] == 0, then we haven't visited the vertex, 
     * if mincut[i] != 0, then the vertex is visited or on the queue. 
     */
    
    /* allocate the queue, it won't ever have more than n entries. */
    q = mxCalloc(sizeof(mbglIndex), n);
    qhead = 0; qtail = 0;
    
    /* add u to the queue */
    mincut[u] = 1;
    q[qtail] = u;
    qtail++;
    
    while (qhead != qtail)
    {
        /* pop the queue */
        v = q[qhead];
        qhead++;
        
        /* look at all neighbors of v. */
        for (k=pia[v]; k < pia[v+1]; k++)
        {        
            /* removed: 2008-04-02: if (cap[k] > 0 && res[k] > 0) */
            if (res[k] > 0)
            {
                /* if res[k] > 0, then the edge is not saturated, so
                 * we follow it find the cut.
                 */
                
                /* get the current neighbor */
                vn = ja[k];
                
                if (mincut[vn] == 0)
                {
                    /* add the vertex to the queue */
                    mincut[vn] = 1;
                    q[qtail] = vn;
                    qtail++;
                }
            }
        }
    }
    
    /* at this point, we've visited all vertices reachable from 
     * unsaturated edges. */
    
    for (k=0; k < n; k++)
    {
        if (mincut[k] == 0)
        {
            mincut[k] = -1;
        }
        
        /* mexPrintf("mincut[%i] = %i\n", k, mincut[k]); */
    }
    
    mxFree (q);
}


/** Check the cut identified against the flow value calculated.
 * 
 * Run through the edges and compute the sum of edges crossing the cut
 * and make sure this is the same as the value of the maximum flow.
 *
 * This function will flag issues with incorrect rounding.
 *
 * @param flow the computed flow value
 * @param pimincut the min cut vector (-1 for side 1, 1 for side 2)
 * @param n, ia, pja, a, the graph structure
 */
void test_cut(int flow, int* pimincut, 
    mbglIndex n, mbglIndex *ia, mbglIndex *pja, double *a)
{
    double cv=0.0; 
    mbglIndex i, j, k;
    for (j=0; j<n; j++) {
        for (k=pja[j]; k<pja[j+1]; k++) {
            i=ia[k];
            /* the source side is 1, and sink side is -1, so we get
             * a source to sink edge if cut[i]>cut[j]. */
            if (pimincut[i]>pimincut[j]) {
                /* mexPrintf("%i -> %i : +%i\n", i+1,j+1,(int)a[k]); */
                cv += a[k];
            }
        }
    }
    if ((int)cv != flow) {
        mexWarnMsgIdAndTxt("max_flow_mex:cutValueNotFlowValue",
          "The rounded (unrounded) value of the minimum cut is %i (%g),"
          "but the value of the max-flow is %i.  These values should be equal",
          (int)cv,cv,flow);
    }
}

/*
 * The mex function runs a max-flow min-cut problem.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mbglIndex i,j,k;
    
    mbglIndex mrows, ncols;
    
    mbglIndex n,nz;
    
    /* sparse matrix */
    mwIndex *A_row, *A_col;
    double *A_val;
    
    /* source/sink */
    mbglIndex u, v;
    
    /* algorithm name */
    char *algname;
    
    /* flow matrix connectivity */
    mbglIndex *pi_flow, *j_flow;
    
    /* capacity and residual structures */
    int *cap, *res;
    
    /* reverse edge map */
    mbglIndex *rev_edge_map;
    
    /* result */
    int flow;
    
    /* output */
    double *pflowval;
    double *pmincut;
    
    double *pri;
    double *prj;
    double *prv;
    
    /* 
     * The current calling pattern is
     * matching_mex(A,verify,initial_match_name,augmenting_path_name)
     */
    
    const mxArray* arg_matrix;
    const mxArray* arg_source;
    const mxArray* arg_sink;
    const mxArray* arg_algname;    
    int required_arguments = 4;
    
    if (nrhs != required_arguments) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the function requires %i arguments, not %i\n", 
            required_arguments, nrhs);
    }
    
    arg_matrix = prhs[0];
    arg_source = prhs[1];
    arg_sink = prhs[2];
    arg_algname = prhs[3];
    
    u = (mbglIndex)load_scalar_arg(arg_source,1);
    v = (mbglIndex)load_scalar_arg(arg_sink,2);
    algname = load_string_arg(arg_algname,3);
    
    /* The first input must be a sparse matrix. */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if (mrows != ncols ||
        !mxIsSparse(prhs[0]) ||
        !mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0])) 
    {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "the matrix must be sparse, square, and double valued");
    }
    
    n = mrows;
    
    /* Get the sparse matrix */
    A_val = mxGetPr(prhs[0]);
    A_row = mxGetIr(prhs[0]);
    A_col = mxGetJc(prhs[0]);
    
    nz = A_col[n];
    
    /* Quick input check */
    if (u > n || u < 1) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "invalid source vertex: %i\n", u);
    } 
    
    if (v > n || v < 1) {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "invalid sink vertex: %i\n", v);
    }
    
    u = u-1;
    v = v-1;
    
    /* build flow connectivity structure */
    build_matrix(n,A_row,A_col,A_val, 
        &pi_flow, &j_flow, &cap, &rev_edge_map);
    
    /* allocate the residual map */
    res = mxCalloc(sizeof(int),pi_flow[n]);
    
    /*i = 0;
    for (k=0; k < pi_flow[n]; k++)
    {
        // get the correct row
        while (k >= pi_flow[i+1]) { ++i; }
        mexPrintf("(%i,%i) (%i,%i)\n", i,j_flow[k],cap[k],res[k]);
    }*/
    
    /* mexPrintf("Calling flow (%i,%i)...\n", u, v); */
    #ifdef _DEBUG
    mexPrintf("max_flow(%s)...", algname);
    #endif 
    if (strcmp(algname,"push_relabel") == 0) {
        push_relabel_max_flow(n,j_flow,pi_flow,
            u,v,cap,res,rev_edge_map,&flow);
    } else if (strcmp(algname, "edmunds_karp") == 0) {
        edmunds_karp_max_flow(n,j_flow,pi_flow,
            u,v,cap,res,rev_edge_map,&flow);
    } else if (strcmp(algname, "kolmogorov") == 0) {
        kolmogorov_max_flow(n,j_flow,pi_flow,
            u,v,cap,res,rev_edge_map,&flow);
    } else {
        mexErrMsgIdAndTxt("matlab_bgl:invalidMexArgument",
            "algname option %s is invalid\n", 
            algname);
    }
    
    #ifdef _DEBUG
    mexPrintf("done!\n");
    #endif 
    
    /*i = 0;
    for (k=0; k < pi_flow[n]; k++)
    {
        // get the correct row
        while (k >= pi_flow[i+1]) { ++i; }
        mexPrintf("(%i,%i) (%i,%i)\n", i,j_flow[k],cap[k],res[k]);
    }*/
    
    if (nlhs >= 1)
    {
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
        pflowval = mxGetPr(plhs[0]);
        pflowval[0] = (double)flow;
    }
    
    if (nlhs >= 2)
    {
        int *pimincut;
        
        plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
        pmincut = mxGetPr(plhs[1]);
        
        pimincut = (int*)pmincut;
        
        build_cut(u, n, pi_flow, j_flow, cap, res, pimincut);
        
        test_cut(flow,pimincut,n,A_row,A_col,A_val);
        
        /* now expand mincut to the full dataset, we need to
         * do this operation backwards because pimincut has integer
         * entries specified and we are expanding them to double.
         */
        expand_int_to_double(pimincut,pmincut,n,0.0);
    }
    
    if (nlhs >= 3)
    {
        plhs[2] = mxCreateDoubleMatrix(nz,1,mxREAL);
        plhs[3] = mxCreateDoubleMatrix(nz,1,mxREAL);
        plhs[4] = mxCreateDoubleMatrix(nz,1,mxREAL);
        
        pri = mxGetPr(plhs[2]);
        prj = mxGetPr(plhs[3]);
        prv = mxGetPr(plhs[4]);
        
        /* j will be our index into the new matrix. */
        j = 0;
        
        for (i=0;i<n;i++)
        {
            for (k=pi_flow[i];k<pi_flow[i+1];k++)
            {
                if (cap[k] != 0)
                {
                    /* since cap[k] != 0, this is a real edge */
                    pri[j] = i+1;
                    prj[j] = j_flow[k]+1;
                    prv[j] = res[k];
                    j++;
                }
            }
        }
        
        if (j != nz)
        {
            mexPrintf("error... j != nz...\n");
        }
    }
    
    #ifdef _DEBUG
    mexPrintf("return\n");
    #endif 
}

