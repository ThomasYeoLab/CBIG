/* 
Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
*/


#include "mex.h"

/* The computational routine */
/*function V_lam=Compute_V_lambda(params,V_same,V_diff,lambda)
neighborsys=double(params.neighborhood);%revise*/
void CBIG_MSHBM_V_lambda_Product(double *neighborsys, double *V_same, double *V_diff, 
    double *lambda, double *V_lam,int N,int K,int M1, int N1)
{
    int m;
    int k;
    int n;
    int j;
    int zj;
    double tmp1;
    double V;
    double V_lam_zj;
    
    /*
    mexPrintf("\nVerts %i ...", N);
    mexPrintf("\nClusters %i ...", K);
    */

    for (m=0; m<N; m++) {
        /*mexPrintf("\nm= %i ...", m);*/
        for (k=0;k<K;k++){
            /*mexPrintf("\nk= %i ...", k);*/
            for (n=0;n<M1;n++){
                /*mexPrintf("\nn= %i ...", n);*/
                j=(int)neighborsys[m*M1+n];
                

                if(j!=0){
                    /*mexPrintf("\nNeighbors %i ...", j);*/
                    
                    for (zj=0;zj<K;zj++){
                        if(zj==k){
                            V=V_same[m*M1+n];
                        }else{
                            V=V_diff[m*M1+n];
                        }

                        tmp1=lambda[zj*N+j-1]*V;
                        if(zj==0){
                            V_lam_zj=tmp1;
                        }else{
                            V_lam_zj=V_lam_zj+tmp1;
                        }
                    }
                     
                }else{
                    /*mexPrintf("\nNeighbors zero ...");*/
                    
                    V_lam_zj=0;
                     
                }
                
                
                
                if(n==0){
                    V_lam[k*N+m]=V_lam_zj;
                }else{
                    V_lam[k*N+m]=V_lam[k*N+m]+V_lam_zj;
                }
                 
               /* V_lam[m*K+k]=lambda[m*K+k];*/    
            }
        }
    }
}

/* The gateway function */
/*(double *neighborsys, double *V_same, double *V_diff, double *lambda, double *V_lam,int m,int n)*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *in_neighborsys;              /* 6xN matrix */
    double *in_V_same;               /* 6xN input matrix */
    double *in_V_diff;               /* 6xN input matrix */
    double *in_lambda;               /* NxK input matrix */
    int ncols;                   /* size of matrix lambda */
    int nrows;
    double *outMatrix;              /* output matrix lambda */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
   /* if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }*/
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    /* check that number of rows in second input argument is 1 */
    /*if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }*/
    
    /* get the value of the scalar input  */
    /*multiplier = mxGetScalar(prhs[0]);*/
    in_neighborsys  = mxGetPr(prhs[0]);
    /* create a pointer to the real data in the input matrix  */
    in_V_same = mxGetPr(prhs[1]);
    in_V_diff = mxGetPr(prhs[2]);
    in_lambda = mxGetPr(prhs[3]);

    /* get dimensions of the input matrix lambda*/
    ncols = mxGetN(prhs[3]);
    nrows = mxGetM(prhs[3]);
    int N1=mxGetN(prhs[0]);
    int M1=mxGetM(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    CBIG_MSHBM_V_lambda_Product(in_neighborsys,in_V_same,in_V_diff,in_lambda,outMatrix,nrows,ncols,M1,N1);
}
