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


void MARS_computeVertexDistSq2Nbors(float *vertexDistSq2Nbors, int maxNeighbors, int numVerts, int *vertexNbors, float *vertices)
{
    int k,j, v_neighbor;
    
    for (k = 0; k < numVerts; k++)
    {
        
        for (j = 0; j < maxNeighbors; j++)
        {
            v_neighbor = vertexNbors[index_2D_array(j, k,maxNeighbors)] - 1; /*convert to C indices*/
            if(v_neighbor != -1)
            {
                /*A real neighbor, compute square distance*/
                vertexDistSq2Nbors[index_2D_array(j, k, maxNeighbors)] = 
                    (vertices[index_2D_array(0,k,3)] - vertices[index_2D_array(0,v_neighbor,3)])*(vertices[index_2D_array(0,k,3)] - vertices[index_2D_array(0,v_neighbor,3)])
                  + (vertices[index_2D_array(1,k,3)] - vertices[index_2D_array(1,v_neighbor,3)])*(vertices[index_2D_array(1,k,3)] - vertices[index_2D_array(1,v_neighbor,3)])
                  + (vertices[index_2D_array(2,k,3)] - vertices[index_2D_array(2,v_neighbor,3)])*(vertices[index_2D_array(2,k,3)] - vertices[index_2D_array(2,v_neighbor,3)]);
                                
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int maxNeighbors, numVerts;
  int *vertexNbors;
  float *vertices;
  int dims[2];
  float *vertexDistSq2Nbors;
  
 /* Check for proper number of arguments. */
  if(nrhs!=4) {
    mexErrMsgTxt("4 inputs required: maxNeighbors (int 32), numVerts (int 32), vertexNbors (int 32 *), vertices (float *)");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /*Get input data*/
  maxNeighbors = * (int *) mxGetData(prhs[0]);
  numVerts = * (int *) mxGetData(prhs[1]);
  vertexNbors = (int *) mxGetData(prhs[2]);
  vertices = (float *) mxGetData(prhs[3]);
  
  /*Allocate Output Data*/
  dims[0] = maxNeighbors;
  dims[1] = numVerts;
  
  /*printf("dims 0 is %d and dims 1 is %d\n", dims[0], dims[1]);*/
  plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
  
  vertexDistSq2Nbors = (float *) mxGetData(plhs[0]);

  /*Fill output Data*/
  MARS_computeVertexDistSq2Nbors(vertexDistSq2Nbors, maxNeighbors, numVerts, vertexNbors, vertices);
  
} 
