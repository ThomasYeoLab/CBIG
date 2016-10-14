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
#include "MARS_linearInterp.h"

int index_2D_array(int row, int col, int num_rows)
{
    return( col*num_rows + row);   
}


void MARS_computeMeshFaceAreas(float *faceAreas, int numFaces, int *faces, float *vertices)
{
    int i, v0, v1, v2;
    
    for (i = 0; i < numFaces; i++)
    {
        v0 = faces[index_2D_array(0, i, 3)] - 1; /*convert to C index*/
        v1 = faces[index_2D_array(1, i, 3)] - 1;
        v2 = faces[index_2D_array(2, i, 3)] - 1;
        
        faceAreas[i] = computeArea(vertices+3*v0, vertices+3*v1, vertices+3*v2);       
    }   
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int numFaces, *faces;
  float *vertices;
  int dims;
  float *faceAreas;
  
 /* Check for proper number of arguments. */
  if(nrhs!=3) {
    mexErrMsgTxt("3 inputs required: numFaces (int 32), faces (int 32 *), vertices (float *)");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /*Get input data*/
  numFaces = * (int *) mxGetData(prhs[0]);
  faces = (int *) mxGetData(prhs[1]);
  vertices = (float *) mxGetData(prhs[2]);
  
  /*Allocate Output Data*/
  dims = numFaces;
  
  plhs[0] = mxCreateNumericArray(1, &dims, mxSINGLE_CLASS, mxREAL);
  
  faceAreas = (float *) mxGetData(plhs[0]);

  /*Fill output Data*/
  MARS_computeMeshFaceAreas(faceAreas, numFaces, faces, vertices);
  
} 
