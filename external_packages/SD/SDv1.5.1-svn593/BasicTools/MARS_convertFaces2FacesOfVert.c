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

/* vertexFaces = MARS_convertFaces2FacesOfVert(faceMatrix, numVertices)

Take in a face matrix (starts count from 1), and the number of vertices
returns a vector (2D matrix) that tells you all the faces a vertex belongs to.
To convert it into a matrix:

num_per_vertex = length(vertexFaces)/size(vertices,1);
vertexFaces = reshape(vertexFaces, size(vertices,1), num_per_vertex);

Now row i of vertexFaces will correspond to the faces vertex i belongs to.
The rows are padded with 0 for vertices with less number of faces than the number of columns of the matrix vertexFaces
*/


/* Assumes faces start with +1 */
mxArray *MARS_convertFaces2FacesOfVert(const mxArray *mxFaces, int numVerts)
{
    int *dimensions;
    int numFaces, *faces;
    int *numFacesPerVertex;
    int maxFaces = 0, row, col;
    mxArray *mxfacesOfVert; 
    int *facesOfVert, vno, total_length;
    
    
    /*First access face matric from mxFaces */
    dimensions = (int *) mxGetDimensions(mxFaces);
    numFaces = dimensions[0];
    faces = (int *) mxGetData(mxFaces);
    
    /* First loop through faces matrix to find the highest number of faces a vertex can be connected to */
    numFacesPerVertex = (int *)calloc(numVerts, sizeof(int));    
    for (row = 0; row < numFaces; row++){
        for(col = 0; col < 3; col++){
            numFacesPerVertex[faces[index_2D_array(row,col,numFaces)] - 1]++; /*-1 because faces start with vertex 1 rather than 0;*/
        }
    }
    
    /*After finishing count for each vertex, now find the maximum number of faces for all the vertices */
    for (row = 0; row < numVerts; row++)
    {
        if(maxFaces < numFacesPerVertex[row])
            maxFaces = numFacesPerVertex[row];
    }
    
    free(numFacesPerVertex);
    
    
    /*Now that we know the maximum number of faces, we can now allocate memory for facesOfVert*/
    total_length = numVerts*maxFaces;
    mxfacesOfVert = mxCreateNumericArray(1, &total_length, mxINT32_CLASS, mxREAL);
    facesOfVert = (int *)mxGetData(mxfacesOfVert);
    
    /*Now we go through the face matrix again, and for each vertex of the face matrix, I assign the corresponding
      face to the vertex*/
    numFacesPerVertex = (int *)calloc(numVerts, sizeof(int)); 
    for (row = 0; row < numFaces; row++){
        for(col = 0; col < 3; col++){
                 
            vno = faces[index_2D_array(row,col,numFaces)]-1; /*subtract 1 because faces start count from vertex 1*/
            
            if(numFacesPerVertex[vno] > maxFaces)
            {
                mexErrMsgTxt("This should not happen!!");
            }
            
            facesOfVert[index_2D_array(vno, numFacesPerVertex[vno], numVerts)] = row + 1; /*+1 because we want to index the faces starting with face 1, not 0*/
            numFacesPerVertex[vno]++;
        }
    }
    /*printf("first entry of facesOfVert is %d\n", facesOfVert[0]);*/
    
    
    free(numFacesPerVertex);
    return(mxfacesOfVert);
}

/*Given row and col of a matrix and the number of rows, return the 1D index of the matrix
  Assume row starts count from 0, col starts count from 0 */
int index_2D_array(int row, int col, int num_rows)
{
    return( col*num_rows + row);   
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int numVerts;

 /* Check for proper number of arguments. */
  if(nrhs!=2) {
    mexErrMsgTxt("Two inputs required. Face matrix and number of vertices, both assumed to be of type int32");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /*plhs[0] = MARS_convertFaces2VertNbors(prhs[0]);*/
  numVerts = * (int *) mxGetData(prhs[1]);
  
  /*printf("num vertices is %d\n", numVerts);*/
  plhs[0] = MARS_convertFaces2FacesOfVert(prhs[0], numVerts);
} 
