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

float triangle_area(float *vertices, int *faces, int face_no, int numFaces, int numVerts)
{
    float vertex1[3];
    float vertex2[3];
    float vertex3[3];
    float u[3];
    float v[3];
    float uv;
    float temp;
    
    vertex1[0] = vertices[index_2D_array(faces[index_2D_array(face_no, 0, numFaces)]-1, 0, numVerts)];
    vertex1[1] = vertices[index_2D_array(faces[index_2D_array(face_no, 0, numFaces)]-1, 1, numVerts)];
    vertex1[2] = vertices[index_2D_array(faces[index_2D_array(face_no, 0, numFaces)]-1, 2, numVerts)];
    
    vertex2[0] = vertices[index_2D_array(faces[index_2D_array(face_no, 1, numFaces)]-1, 0, numVerts)];
    vertex2[1] = vertices[index_2D_array(faces[index_2D_array(face_no, 1, numFaces)]-1, 1, numVerts)];
    vertex2[2] = vertices[index_2D_array(faces[index_2D_array(face_no, 1, numFaces)]-1, 2, numVerts)];
    
    vertex3[0] = vertices[index_2D_array(faces[index_2D_array(face_no, 2, numFaces)]-1, 0, numVerts)];
    vertex3[1] = vertices[index_2D_array(faces[index_2D_array(face_no, 2, numFaces)]-1, 1, numVerts)];
    vertex3[2] = vertices[index_2D_array(faces[index_2D_array(face_no, 2, numFaces)]-1, 2, numVerts)];
    
    u[0] = vertex2[0] - vertex1[0];
    u[1] = vertex2[1] - vertex1[1];
    u[2] = vertex2[2] - vertex1[2];
    
    v[0] = vertex3[0] - vertex1[0];
    v[1] = vertex3[1] - vertex1[1];
    v[2] = vertex3[2] - vertex1[2];
    
    /* 0.5 * | u x v| = 0.5 * |u||v|sin(theta) = 0.5 * sqrt(|u|^2|v|^2 - (u.v)^2) */
    temp = (u[0]*u[0]+u[1]*u[1]+u[2]*u[2])*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) - (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])*(u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
    
    if(temp < 0)
    {  
        temp = 0.0;
    }
    return(0.5 * sqrt( temp ));
    
}


/* Assumes faces start with +1 */
float MARS_calculateSurfaceArea(const mxArray *mxVerts, const mxArray *mxFaces)
{
    int *dimensions;
    int numFaces, numVerts;
    int *faces;
    float *vertices;
    int i; 
    float surfaceArea = 0.0;
    
    dimensions = (int *) mxGetDimensions(mxFaces);
    numFaces = dimensions[0];
    faces = (int *) mxGetData(mxFaces);
    
    dimensions = (int *) mxGetDimensions(mxVerts);
    numVerts = dimensions[0];
    vertices = (float *) mxGetData(mxVerts);
    
    for (i = 0; i < numFaces; i++)
    {
        surfaceArea += triangle_area(vertices, faces, i, numFaces, numVerts);
    }
    
    return(surfaceArea);
}

/* Assume row starts count from 0, col starts count from 0 */
int index_2D_array(int row, int col, int num_rows)
{
    return( col*num_rows + row);   
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize length = 1;
  float *surfaceArea;
  
 /* Check for proper number of arguments. */
  if(nrhs!=2) {
    mexErrMsgTxt("Two inputs required. Face matrix and number of vertices, both assumed to be of type int32");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  mexErrMsgTxt("There is a bug in this code. Use sum(MARS_computeMeshFaceAreas) instead");
  
  plhs[0] = mxCreateNumericArray(1, &length, mxSINGLE_CLASS, mxREAL);
  surfaceArea = (float *) mxGetData(plhs[0]);
  
  *surfaceArea = MARS_calculateSurfaceArea(prhs[0], prhs[1]);
} 
