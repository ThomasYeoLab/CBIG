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

	File:	MARS_linearInterpolateAux.c
	A mex file for Matlab (c)
	The calling syntax is:
	
    [val, grad, NViF] = MARS_linearInterpolateAuxWGrad(points, mesh.vertices, mesh.faces, mesh.vertexNbors, 
                                      mesh.vertexFaces, seedVertices, data)
    
	 points = 3 x numOfPoints (single precision)
     vertices = 3 x numOfVertices (single)
     faces = 3 x numOfTriangles (int32)
 *   vertexNbors = maxNumOfNbors x numOfVertices (int32)
 *   vertexFaces = maxNumOfFaces x numOfVertices (int32)
 *   data = d x numOfPoints (single);
 *   seedVertices = 1 x numOfVertices (int32)   
 *
 *   [val] = d x numOfPoints (single) 
 *   [grad] = d x 3 x numOfPoints (single)
 *   [NViF] = 1 x numOfPoints (int32)
	 Date:	October 2006

	Copyright (c) 2006 by Mert Rory Sabuncu
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
 * The main routine.  
 */
void
mexFunction(
			int nlhs, 
			mxArray *plhs[], 
			int nrhs, 
			const mxArray *prhs[])
{

    /*
     points = 3 x numOfPoints (single precision)
     vertices = 3 x numOfVertices (single)
     faces = 3 x numOfTriangles (int32)
 *   vertexNbors = maxNumOfNbors x numOfVertices (int32)
 *   vertexFaces = maxNumOfFaces x numOfVertices (int32)
 *   data = d x numOfVertices (single);
 *   seedVertices = 1 x numOfVertices (int32)  
    */
    
    const float *	points   		=	(float*) mxGetPr(prhs[0]);
	const float *   vertices		=	(float*) mxGetPr(prhs[1]);
    const int   *   faces           =   (int*)  mxGetPr(prhs[2]);
	const int	*	vertexNbors	    =   (int*)  mxGetPr(prhs[3]);
    const int   *   vertexFaces     =   (int*)  mxGetPr(prhs[4]);
    const int   *   seedVertices    =   (int*)  mxGetPr(prhs[5]);
    const float *   data            = (float*)  mxGetPr(prhs[6]);
    
    const int       numOfNbors      =   mxGetM(prhs[3]);
    const int       numOfVertFaces  =   mxGetM(prhs[4]);
    const int       numOfVertices   =   mxGetN(prhs[1]);
    const int       numOfFaces      =   mxGetN(prhs[2]);
    const int       numOfPoints     =   mxGetN(prhs[0]);
    const int       data_dim        =   mxGetM(prhs[6]);
    
    int out_dims[2], p, grad_dims[3];
    int index_dims[2];
    int * FaceInds_MatlabStyle, * NViF_MatlabStyle;
    int cur_face_ind;
    
    const float * v0;
    const float * v1;
    const float * v2;
    const float * d0;
    const float * d1;
    const float * d2;
    const float * cur_point;
    float * new_data;
    float * grad;
    
    out_dims[0] = data_dim;
    out_dims[1] = numOfPoints;
    grad_dims[0] = data_dim;
    grad_dims[1] = 3;
    grad_dims[2] = numOfPoints;
    index_dims[0] = 1;
    index_dims[1] = numOfPoints;
    
    #ifdef DEBUG
    
    printf("NumOfPoints : %i, Num of Faces: %i, Num of Vertices: %i, Num of Max Nbors: %i, Num of Max Vert Faces: %i\n",
            numOfPoints, numOfFaces, numOfVertices, numOfNbors, numOfVertFaces);
    
    #endif
    
    if (mxGetM(prhs[0]) != 3)	mexErrMsgTxt("Point should be three dimensional!");
    if (mxGetN(prhs[1]) != mxGetN(prhs[3]) 
        && mxGetN(prhs[1]) != mxGetN(prhs[4])) 
        mexErrMsgTxt("Inconsistent number of vertices!");
    if (mxGetN(prhs[0]) != mxGetN(prhs[5]))
        mexErrMsgTxt("Numof points and number of seedVertices should be the same!");
    if (mxGetN(prhs[6]) != numOfVertices) mxErrMsgTxt("Data not the right size!");
    
    /*
     *
     *First find faces
     */
    
    /*if(!(FaceInds_MatlabStyle = (int*) calloc(numOfPoints, sizeof(int))))
         mxErrMsgTxt("Memory allocation error!!!");*/
    plhs[3] = mxCreateNumericArray(2, index_dims, mxINT32_CLASS, mxREAL);
    FaceInds_MatlabStyle = (int *)mxGetData(plhs[3]);    
     
    plhs[2] = mxCreateNumericArray(2, index_dims, mxINT32_CLASS, mxREAL);
    NViF_MatlabStyle = (int*) mxGetData(plhs[2]);
    
    
    MARS_findFaces(points, numOfPoints, vertices, numOfVertices, faces, numOfFaces,
                vertexNbors, numOfNbors, vertexFaces, numOfVertFaces, seedVertices, FaceInds_MatlabStyle, NViF_MatlabStyle);
    
    plhs[0] = mxCreateNumericArray(2, out_dims, mxSINGLE_CLASS, mxREAL);    
    new_data = (float *) mxGetData(plhs[0]);    
     
    plhs[1] = mxCreateNumericArray(3, grad_dims, mxSINGLE_CLASS, mxREAL);    
    grad = (float *) mxGetData(plhs[1]);
    
    
    for (p = 0; p < numOfPoints; p++)
    {
        /*
         *Convert to C-style
         */       
        
        cur_point = &points[3 * p];
        cur_face_ind = FaceInds_MatlabStyle[p] - 1;
        v0 = &vertices[3*(faces[3*cur_face_ind] - 1)];
        v1 = &vertices[3*(faces[3*cur_face_ind + 1] - 1)];
        v2 = &vertices[3*(faces[3*cur_face_ind + 2] - 1)];
        
        d0 = &data[data_dim *(faces[3*cur_face_ind] - 1)];
        d1 = &data[data_dim *(faces[3*cur_face_ind + 1] - 1)];
        d2 = &data[data_dim *(faces[3*cur_face_ind + 2] - 1)];
        
        #ifdef DEBUG
        printf("Interpolating at point: %f %f %f\n", cur_point[0], cur_point[1], cur_point[2]);
        printf("v0: %f %f %f\n", v0[0], v0[1], v0[2]);
        printf("v1: %f %f %f\n", v1[0], v1[1], v1[2]);
        printf("v2: %f %f %f\n", v2[0], v2[1], v2[2]);
        printf("Data:\n");
        printf("d0: %f %f %f\n", d0[0], d0[1], d0[2]);
        printf("d1: %f %f %f\n", d1[0], d1[1], d1[2]);
        printf("d2: %f %f %f\n", d2[0], d2[1], d2[2]);
        #endif
        linearInterpWGrad(cur_point, v0, v1, v2, d0, d1, d2, data_dim, &(new_data[data_dim * p]), &(grad[data_dim * 3 * p]));   
    } 
    
    /*free(FaceInds_MatlabStyle);*/
};
