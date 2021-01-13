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
 * at each mesh vertex, assuming a displacement vector at each mesh vertex. 
 * It asssumes barycentric interpolation followed by normalization to the sphere. 
 * This is done by a weighted average of the gradient of all the 
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
    const int       geodesic        = * (int *) mxGetData(prhs[5]);
    
    const int       numOfVertFaces  =   mxGetM(prhs[2]);
    const int       numOfVertices   =   mxGetN(prhs[0]);
    const int       data_dim        =   mxGetM(prhs[4]);
    
    
    const float radius_sq = dotVectors(vertices, vertices);    
    float Dx[9], Dx_sq[9], temp_grad2[9];
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
    
    if(data_dim != 3)
        mexErrMsgTxt("Assume data dim is 3!!"); 
       
    
    
    plhs[0] = mxCreateNumericArray(3, grad_dims, mxSINGLE_CLASS, mxREAL);    
    grad = (float *) mxGetData(plhs[0]);
    temp_grad = (float *) calloc(data_dim*3, sizeof(float));
    new_data = (float *) calloc(data_dim, sizeof(float));
    
    /*initialize to 0.*/
    for (i = 0; i < data_dim*3*numOfVertices; i++)
       grad[i] = 0;
    
    
    /*printf("%f\n", radius_sq);*/
    for (p = 0; p < numOfVertices; p++) /*shoudl be numOfVertices*/
    {
        /*
         *Convert to C-style
         */       
        
        cur_vertex = &vertices[3 * p];
        total_area = 0;
        
        
        if(geodesic != 1)
        {
        /*
         * Compute 1/r^2 * Dx^2
         */
            
            Dx[0] = 0;
            Dx[1] = cur_vertex[2];
            Dx[2] = -cur_vertex[1];
            
            Dx[3] = -cur_vertex[2];
            Dx[4] = 0;
            Dx[5] = cur_vertex[0];
            
            Dx[6] = cur_vertex[1];
            Dx[7] = -cur_vertex[0];
            Dx[8] = 0;
            
            MultiplyMatrix(Dx, Dx, Dx_sq, 3, 3, 3);
            for (i = 0; i < 9; i++)
                Dx_sq[i] = Dx_sq[i]/radius_sq;
        }
        else
        {
            Dx_sq[0] = 1;
            Dx_sq[1] = 0;
            Dx_sq[2] = 0;
            
            Dx_sq[3] = 0;
            Dx_sq[4] = 1;
            Dx_sq[5] = 0;
            
            Dx_sq[6] = 0;
            Dx_sq[7] = 0;
            Dx_sq[8] = 1;
        }
        
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
         
                 barycentricInterpWGradWONormalization(cur_vertex, v0, v1, v2, d0, d1, d2, data_dim, new_data, temp_grad);   
                 
                 for (i = 0; i < data_dim*3; i++)
                 {
                     grad[p*data_dim*3 + i] += temp_grad[i] * faceAreas[fno_c_style];
                 }
            }
        }
        for (i = 0; i < data_dim*3; i++)
        {
            temp_grad2[i] = grad[p*data_dim*3 + i]/total_area;
        }
        
/*         if(p == 0){
         printf("========\n");       
         printf("%f %f %f\n", temp_grad2[0], temp_grad2[3], temp_grad2[6]);
         printf("%f %f %f\n", temp_grad2[1], temp_grad2[4], temp_grad2[7]);
         printf("%f %f %f\n", temp_grad2[2], temp_grad2[5], temp_grad2[8]);
         printf("========\n");
         }
         
         if(p == 0){
         printf("========\n");       
         printf("%f %f %f\n", Dx_sq[0], Dx_sq[3], Dx_sq[6]);
         printf("%f %f %f\n", Dx_sq[1], Dx_sq[4], Dx_sq[7]);
         printf("%f %f %f\n", Dx_sq[2], Dx_sq[5], Dx_sq[8]);
         printf("========\n");
         }*/
        
        MultiplyMatrix(Dx_sq, temp_grad2, &grad[p*data_dim*3], 3, 3, 3);       
    } 
 /*   
     printf("========\n");       
     printf("%f %f %f\n", grad[0], grad[3], grad[6]);
     printf("%f %f %f\n", grad[1], grad[4], grad[7]);
     printf("%f %f %f\n", grad[2], grad[5], grad[8]);
     printf("========\n");*/
    free(temp_grad);
    free(new_data);
};
