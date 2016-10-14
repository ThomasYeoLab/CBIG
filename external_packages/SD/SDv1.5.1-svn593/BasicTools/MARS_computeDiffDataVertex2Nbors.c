/* Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions */

/*=========================================================================

  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the names of the copyright holders nor the names of future
      contributors may be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    

=========================================================================*/
#include "mex.h"
#include "math.h"
#include "matrix.h"

int index_2D_array(int row, int col, int num_rows)
{
    return( col*num_rows + row);   
}

/* Assume all dimensions start count from  0 */
int index_3D_array(int dim1, int dim2, int dim3, int num_dim1, int num_dim2)
{
    return(dim3*(num_dim1*num_dim2) + dim2*num_dim1 + dim1);
}

void MARS_computeDiffDataVertex2Nbors(float *diffDataVertex2Nbors, int maxNeighbors, int numVerts, int numData, int *vertexNbors, float *data)
{
    int k,j, d, v_neighbor;
    
    for (k = 0; k < numVerts; k++) /*for each vertex*/
    {     
        for (j = 0; j < maxNeighbors; j++) /*for each possible neighbor*/
        {
            v_neighbor = vertexNbors[index_2D_array(j, k,maxNeighbors)] - 1; /*convert to C indices*/
            if(v_neighbor != -1) /*real neighbor*/
            {
                for (d = 0; d < numData; d++) /* for each data dimension*/
                {
                    diffDataVertex2Nbors[index_3D_array(d, j, k, numData, maxNeighbors)] = data[index_2D_array(d, k, numData)] - data[index_2D_array(d, v_neighbor,numData)]; 
                }
            }
        }
    }
    
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int maxNeighbors, numVerts, numData;
  int *vertexNbors;
  float *data;
  int dims[3];
  float *diffDataVertex2Nbors;
  
 /* Check for proper number of arguments. */
  if(nrhs!=5) {
    mexErrMsgTxt("5 inputs required: maxNeighbors (int 32), numVerts (int 32), numData (int 32), vertexNbors (int 32 *), data (float *)");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /*Get input data*/
  maxNeighbors = * (int *) mxGetData(prhs[0]);
  numVerts = * (int *) mxGetData(prhs[1]);
  numData = * (int *) mxGetData(prhs[2]);
  vertexNbors = (int *) mxGetData(prhs[3]);
  data = (float *) mxGetData(prhs[4]);
  
  /*Allocate Output Data*/
  dims[0] = numData;
  dims[1] = maxNeighbors;
  dims[2] = numVerts;
  
  plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL); /*numData x maxNeighbors x numVerts*/
  diffDataVertex2Nbors = (float *) mxGetData(plhs[0]);

  /*Fill output Data*/
  MARS_computeDiffDataVertex2Nbors(diffDataVertex2Nbors, maxNeighbors, numVerts, numData, vertexNbors, data);
  
} 
