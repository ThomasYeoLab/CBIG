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
#include <memory.h>

/* Assume row starts count from 0, col starts count from 0 */
int index_2D_array(int row, int col, int num_rows)
{
    return( col*num_rows + row);   
}

void Find1NeighborhoodGivenRadiusAux(int maxNeighbors, int numVerts, int *vertexNbors, int *indicator_vec1, int *indicator_vec2)
{
    int i, j, vno;
    
    for(i = 0; i < numVerts; i++) /*for each vertex*/
    {
        if(indicator_vec1[i] == 1) /*if vertex is in neighborhood*/
        {
            for(j = 0; j < maxNeighbors; j++) /*check its neighbors*/
            {
                vno = vertexNbors[index_2D_array(j, i, maxNeighbors)]; /*if real neighbor*/
                if(vno != 0)
                    indicator_vec2[vno-1] = 1; /*add neighbor to the list. subtract 1 for C-index*/
            }
        }
    }
    
    memcpy(indicator_vec1, indicator_vec2, numVerts*sizeof(int));
}

void FindNeighborhoodGivenRadiusAux(mxArray *plhs[], int maxNeighbors, int numVerts, int vertex, int radius, int *vertexNbors)
{
    int *indicator_vec1, *indicator_vec2;
    int i, numVertInNeighborhood = 0;
    int dims[2];
    int *vertexList;
    
    indicator_vec1 = (int *) calloc(numVerts, sizeof(int));
    indicator_vec2 = (int *) calloc(numVerts, sizeof(int));
    
    /*Initialization.*/
    indicator_vec1[vertex-1] = 1; /*-1 to convert to C index.*/
    indicator_vec2[vertex-1] = 1;
    
    /*For each radius, expand it once*/
    for (i = 0; i < radius; i++)
        Find1NeighborhoodGivenRadiusAux(maxNeighbors, numVerts, vertexNbors, indicator_vec1, indicator_vec2);
    
    /*Count number of vertices in neighborhood*/
    for (i = 0; i < numVerts; i++)
    {
        if(indicator_vec1[i] == 1)
            numVertInNeighborhood++;
    }
    
    /*allocate space for output*/
    dims[0] = 1;
    dims[1] = numVertInNeighborhood;
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    vertexList = (int *) mxGetData(plhs[0]);
    
    
    numVertInNeighborhood = 0;
    for (i = 0; i < numVerts; i++) /*for each vertex*/
    {
        if(indicator_vec1[i] == 1) /*if vertex is in neighborhood*/
        {
            vertexList[numVertInNeighborhood] = i+1; /*convert to matlab index and add to neighborhood list*/
            numVertInNeighborhood++;
        }
    }
    free(indicator_vec1);
    free(indicator_vec2);
    
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int numVerts;
  int maxNeighbors;
  int vertex;
  int radius;
  int *vertexNbors;
  
 /* Check for proper number of arguments. */
  if(nrhs!=5) {
    mexErrMsgTxt("Five inputs required. maxNeighbors (int32), numVerts (int32), vertex (int32), radius (int32), vertexNbors (int *)");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }

  /*Get input data*/
  maxNeighbors = * (int *) mxGetData(prhs[0]);
  numVerts = * (int *) mxGetData(prhs[1]);
  vertex = * (int *) mxGetData(prhs[2]);
  radius = * (int *) mxGetData(prhs[3]);
  vertexNbors =  (int *) mxGetData(prhs[4]);
  
  /*Fill output Data*/
  FindNeighborhoodGivenRadiusAux(plhs, maxNeighbors, numVerts, vertex, radius, vertexNbors);
} 




