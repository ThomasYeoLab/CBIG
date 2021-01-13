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
/***********************************************************************

	File:	MARS_linearInterpolateVertexAuxWGrad.c
	A mex file for Matlab (c)
	The calling syntax is:
	
    [grad] = MARS_linearInterpolateVertexAuxWGrad(mesh.vertices, mesh.faces, mesh.vertexNbors, 
                                      mesh.vertexFaces, mesh.faceAreas, data)
    
     vertices = 3 x numOfVertices (single)
     faces = 3 x numOfTriangles (int32)
 *   vertexNbors = maxNumOfNbors x numOfVertices (int32)
 *   vertexFaces = maxNumOfFaces x numOfVertices (int32)
 *   faceAreas = 1 x numOfFaces (single)
     data = d x numOfPoints (single);
 * 
 *   [grad] = d x 3 x numOfPoints (single)
 
 Date:	Feb 2008

	Copyright (c) 2008 by Thomas Yeo
************************************************************************/

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#define MAX_DISTANCE_SQ 100000
#define MIN_TOL     0.00005
/*#define DEBUG*/
#include "MARS_vec3D.h"
#include "MARS_linearInterp.h"
#include "MARS_findFaces.h"


/*
 * The main routine. This function finds the gradient of a mesh
 * at each mesh vertex, assuming scalar values at each vertex. 
 * It assumes barycentric interpolation and weights the gradient of all the 
 * faces connected to a vertex.
 */
void
mexFunction(
			int nlhs, 
			mxArray *plhs[], 
			int nrhs, 
			const mxArray *prhs[])
{

    /*
     * vertices = 3 x numOfVertices (single)
     * faces = 3 x numOfTriangles (int32)
     * vertexNbors = maxNumOfNbors x numOfVertices (int32)
     * vertexFaces = maxNumOfFaces x numOfVertices (int32)
     * faceAreas = 1 x numOfFaces (single)
     * data = d x numOfPoints (single);
     * 
     * [grad] = d x 3 x numOfPoints (single)
     *
     */
    
    const float *   vertices		=  (float*) mxGetPr(prhs[0]);
    const int   *   faces           =   (int*)  mxGetPr(prhs[1]);
	const int   *   vertexFaces     =   (int*)  mxGetPr(prhs[2]);
    const float *   faceAreas       = (float*)  mxGetPr(prhs[3]);
    const float *   data            = (float*)  mxGetPr(prhs[4]);
    
    const int       numOfVertFaces  =   mxGetM(prhs[2]);
    const int       numOfVertices   =   mxGetN(prhs[0]);
    const int       data_dim        =   mxGetM(prhs[4]);
    
    mwSize grad_dims[3];
    int p, f, fno_c_style;
    int i;
    
    const float * v0;
    const float * v1;
    const float * v2;
    const float * d0;
    const float * d1;
    const float * d2;
    const float * cur_vertex;
    float * new_data;
    float * grad;
    float * temp_grad;
    float total_area;
    
    grad_dims[0] = data_dim;
    grad_dims[1] = 3;
    grad_dims[2] = numOfVertices;
    
    
    plhs[0] = mxCreateNumericArray(3, grad_dims, mxSINGLE_CLASS, mxREAL);    
    grad = (float *) mxGetData(plhs[0]);
    temp_grad = (float *) calloc(data_dim*3, sizeof(float));
    new_data = (float *) calloc(data_dim, sizeof(float));
    
    /*initialize to 0.*/
    for (i = 0; i < data_dim*3*numOfVertices; i++)
       grad[i] = 0;
    
    
    
    
    for (p = 0; p < numOfVertices; p++) /*shoudl be numOfVertices*/
    {
        /*
         *Convert to C-style
         */       
        
        cur_vertex = &vertices[3 * p];
        total_area = 0;
        
        for (f = 0; f < numOfVertFaces; f++) /*numOfVertFaces*/
        {
             fno_c_style = vertexFaces[p*numOfVertFaces + f] - 1;
             
             if(fno_c_style != -1)
             {
                 total_area += faceAreas[fno_c_style];
      
                 v0 = &vertices[3*(faces[3*fno_c_style] - 1)];
                 v1 = &vertices[3*(faces[3*fno_c_style + 1] - 1)];
                 v2 = &vertices[3*(faces[3*fno_c_style + 2] - 1)];
             
                 d0 = &data[data_dim *(faces[3*fno_c_style] - 1)];
                 d1 = &data[data_dim *(faces[3*fno_c_style + 1] - 1)];
                 d2 = &data[data_dim *(faces[3*fno_c_style + 2] - 1)];
         

                 
                 barycentricInterpWGrad(cur_vertex, v0, v1, v2, d0, d1, d2, data_dim, new_data, temp_grad);   
                 
                 for (i = 0; i < data_dim*3; i++)
                 {
                     grad[p*data_dim*3 + i] += temp_grad[i] * faceAreas[fno_c_style];
                     
                 }
                 
            }
        }
        for (i = 0; i < data_dim*3; i++)
        {
            grad[p*data_dim*3 + i] = grad[p*data_dim*3 + i]/total_area;
        }
        
    } 
    
    free(temp_grad);
    free(new_data);
};
