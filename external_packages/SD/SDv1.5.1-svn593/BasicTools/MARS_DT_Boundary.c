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
#include "min_heap.h"

#define BIG_NO 100000000.0

int index_2D_array(int row, int col, int num_rows)
{
    return( col*num_rows + row);   
}

void MARS_DT_Boundary(int *sv, int numVerts, int maxNeighbors, int *vertNbors, double *costNbors, double *final_cost)
{
    MIN_HEAP *MH;
    int *prev_neighbor, neighbor;
    int i, result, heap_size;
    double key; int id; void *data; double queryKey;
    double cost;
    
    prev_neighbor = (int *) calloc(numVerts, sizeof(int)); /*Note that this is initialized to 0*/
    MH = Min_HeapAllocate(numVerts, numVerts);
    
    for(i=0; i<numVerts; i++)
    {
        result = Min_HeapInsert(MH, BIG_NO, NULL, (int) i); /*C indices*/
        if(result == ERROR)
            mexErrMsgTxt("MARS_DT_Boundary: insert outside range!!!"); 
        
        /*Initialize sources (boundaries) to be zeros!*/
        if(sv[i] == 1)
            result = Min_HeapEditKeyIndexID(MH, (int) i, 0.0); /*C indices*/
            
        if(result == ERROR)
            mexErrMsgTxt("MARS_DT_Boundary: edit outside range!!!");    
    }
    
    while(1)
    {
        heap_size = Min_HeapGetCurrSize(MH);
        if(heap_size == 0)
        {
            break; /*Completed distance transform!*/
        }
        
        result = Min_HeapExtract(MH, &key, &data, &id);
        if(result == ERROR)
            mexErrMsgTxt("MARS_DT_Boundary: extract fail!!!");
        
        final_cost[id] = key;
        for(i = 0; i < maxNeighbors; i++)
        {
            neighbor = vertNbors[index_2D_array(i, id, maxNeighbors)]; /*id is already C-index*/
            if(neighbor != 0) /*real neighbor*/
            {
                if(neighbor > numVerts || neighbor < 1)
                        mexErrMsgTxt("MARS_DT_Boundary: neighbor outside range!!!");
                
                neighbor--; /*convert to C index*/
                if(Min_HeapIdIsInHeap(MH, neighbor)) /*Neighbor has not been taken out of heap yet*/
                {
                    result = Min_HeapQueryKeyIndexID(MH, neighbor, &queryKey);
                    if(result == ERROR)
                          mexErrMsgTxt("MARS_DT_Boundary: query fail!!!");
                    
                    cost = costNbors[index_2D_array(i, id, maxNeighbors)]; /*note that id is already in C index*/
                    if(queryKey > key + cost)
                        {
                            result = Min_HeapEditKeyIndexID(MH, neighbor, key+cost);
                            if(result == ERROR)
                                mexErrMsgTxt("dijkstra_src2dest: edit fail!!!");
                            prev_neighbor[neighbor] = id; /*neighbor already at C index*/
                        }                   
                }               
                
            }
        }       
    }
    free(prev_neighbor);
    Min_HeapFree(MH);
    
    return;
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{      
   double *final_cost;
   int *sv;
   int *vertNbors; 
   double *costNbors;
   int numVerts, maxNeighbors;
   int dims[1];

   /* Assume inputs are all in matlab indices
    Check for proper number of arguments. */
   if(nrhs!=5) {
      mexErrMsgTxt("Six inputs required.");
   } else if(nlhs>2) {
      mexErrMsgTxt("Too many output arguments");
   }

   sv = (int *) mxGetData(prhs[0]);
   numVerts = * (int *) mxGetData(prhs[1]); 
   maxNeighbors = * (int *) mxGetData(prhs[2]); 
   vertNbors = (int *) mxGetData(prhs[3]);
   costNbors = (double *) mxGetData(prhs[4]); /*unfortunately min_heap assumes double.*/
   
   dims[0] = numVerts;
   plhs[0] = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxREAL);
   final_cost = (double *) mxGetPr(plhs[0]);
   
   MARS_DT_Boundary(sv, numVerts, maxNeighbors, vertNbors, costNbors, final_cost);    
} 


