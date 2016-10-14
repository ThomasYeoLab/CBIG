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

/* vertexFaces = MARS_convertFaces2VertNbors(faceMatrix, numVertices)

Take in a face matrix (starts count from 1), and the number of vertices
returns a vector (2D matrix) that tells you all the neighbors of a vertex.
To convert it into a matrix:

num_per_vertex = length(vertexNbors)/size(vertices,1);
vertexNbors = reshape(vertexNbors, size(vertices,1), num_per_vertex);

Now row i of vertexNbors will correspond to the neighbors of vertex i
The rows are padded with 0 for vertices with less neighbors than the number of columns of the matrix vertexFaces
*/

/* Adding neighbors. Assume vno1 and vno2 starts from 0.
   Given a pair of vertices which are neighbors. Add them as neighbors to each other if necessary (avoid duplication)*/
void AddNeighbors(int vno1, int vno2, int *vertNbors, int *numNborsPerVertex, int numVerts, int maxFaces)
{
    int found=0;
    int i;
    
    found = 0;  /*First check whether the vertices are already connected. */
    for(i = 0; i < numNborsPerVertex[vno1]; i++){
        
        if(vno2 == (vertNbors[index_2D_array(vno1, i, numVerts)]-1)){ /*-1 because we save vertices as starting from +1, but vno1 and vno2 starts at 0.*/
            found = 1; 
            break;
        }
    }
    if(found==0){
        /*Not added before, can add it to both.*/
        vertNbors[index_2D_array(vno1, numNborsPerVertex[vno1], numVerts)] = vno2+1; /*+1 so that we start with vertex 1*/
        numNborsPerVertex[vno1]++;
        vertNbors[index_2D_array(vno2, numNborsPerVertex[vno2], numVerts)] = vno1+1; /*+1 so that we start with vertex 1*/
        numNborsPerVertex[vno2]++;
    }       
}


mxArray *MARS_convertFaces2VertNbors(const mxArray *mxFaces, int numVerts)
{
    int *dimensions;
    int numFaces, *faces;
    int *numFacesPerVertex, row, col, maxFaces = 0;
    
    int face_no, vno1, vno2, vno3, i;
    mxArray *mxVertNbors; int *vertNbors, total_length;
    int *numNborsPerVertex;
    
    /*First access the face matrix*/
    dimensions = (int *) mxGetDimensions(mxFaces);
    /*printf("faces: dimension 1 is %d and dimension 2 is %d\n", dimensions[0], dimensions[1]);*/
    
    numFaces = dimensions[0];
    faces = (int *) mxGetData(mxFaces);
    
    
    /* First loop through faces matrix to find the highest number of faces a vertex can be connected to*/
    numFacesPerVertex = (int *)calloc(numVerts, sizeof(int));    
    for (row = 0; row < numFaces; row++){
        for(col = 0; col < 3; col++){
            numFacesPerVertex[faces[index_2D_array(row,col,numFaces)] - 1]++; /*-1 because faces start with 1 rather than 0;*/
        }
    }
    
    /*After finishing count for each vertex, now find the maximum number of faces for all the vertices*/
    for (row = 0; row < numVerts; row++)
    {
        if(maxFaces < numFacesPerVertex[row])
            maxFaces = numFacesPerVertex[row];
    }
    
    /*printf("Max faces is %d\n", maxFaces);*/
    free(numFacesPerVertex);
    
    /*Now that we know the max number of faces, we can allocate memory. Note the number of neighbors a vertex has = the number of faces a vertex is involved in (for a spherical mesh)*/
    total_length = numVerts*maxFaces;
    mxVertNbors = mxCreateNumericArray(1, &total_length, mxINT32_CLASS, mxREAL);
    vertNbors = (int *) mxGetData(mxVertNbors); 
    numNborsPerVertex = (int *)calloc(numVerts, sizeof(int));
    
    /*Now for each face, we extract all the vertices, and we add the vertices of each other if they have not been already added (avoid duplication)*/
    for (face_no = 0; face_no < numFaces; face_no++){
        
        vno1 = faces[index_2D_array(face_no,0,numFaces)]-1;
        vno2 = faces[index_2D_array(face_no,1,numFaces)]-1;
        vno3 = faces[index_2D_array(face_no,2,numFaces)]-1;
        
        AddNeighbors(vno1, vno2, vertNbors, numNborsPerVertex, numVerts, maxFaces);
        AddNeighbors(vno1, vno3, vertNbors, numNborsPerVertex, numVerts, maxFaces);
        AddNeighbors(vno2, vno3, vertNbors, numNborsPerVertex, numVerts, maxFaces);        
    }
    
    free(numNborsPerVertex);
    return(mxVertNbors);
}



/* Assume row starts count from 0, col starts count from 0 */
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
  plhs[0] = MARS_convertFaces2VertNbors(prhs[0], numVerts);
} 
