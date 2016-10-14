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

	File:	MARS_findFaceAux.c
	A mex file for Matlab (c)
	The calling syntax is:
			
	 [FaceInds] = MARS_findFaceAux(points, mesh.vertices, mesh.faces, mesh.vertexNbors, mesh.vertexFaces, seedVertices)
	
	 point = 3 x numOfPoints (single precision)
     mesh.vertices = 3 x numOfVertices (single) 
 *   mesh.faces    = 3 x numOfFaces    (int32)
	 mesh.vertexNbors = N x numOfVertices (int32)
 *   mesh.vertexFaces = M x numOfVertices (int32)
     seedVertices = 1 x numOfPoints (int32)
     FaceInd = 1 x numOFPoints (int32) (NOTE: Matlab-style index!!!)
 *
 *   nearestVertInFace = 1 x numOfPoints (int32)
      
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
#include "MARS_vec3D.h"

/*#define DEBUG2
  #define DEBUG
  #define DEBUG3*/

#define MAX_DISTANCE_SQ 100000
#define MIN_TOL     0.00001

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
	/*Declarations*/
    /*Inputs:
     points = 3 x numOFPoints (single precision)
     mesh.vertices = 3 x numOfVertices (single) 
 *   mesh.faces    = 3 x numOfFaces    (int32)
	 mesh.vertexNbors = N x numOfVertices (int32)
 *   mesh.vertexFaces = M x numOfVertices (int32)
     seedVertex = scalar (int32)
    */
	
	const float *	points   		=	(float*) mxGetPr(prhs[0]);
	const float *   vertices		=	(float*) mxGetPr(prhs[1]);
    const int   *   faces           =   (int*)  mxGetPr(prhs[2]);
	const int	*	vertexNbors	    =   (int*)  mxGetPr(prhs[3]);
    const int   *   vertexFaces     =   (int*)  mxGetPr(prhs[4]);
    const int   *   seedVertices    =   (int*)  mxGetPr(prhs[5]);
    
    const int       numOfNbors      =   mxGetM(prhs[3]);
    const int       numOfVertFaces  =   mxGetM(prhs[4]);
    const int       numOfVertices   =   mxGetN(prhs[1]);
    const int       numOfFaces      =   mxGetN(prhs[2]);
    const int       numOfPoints     =   mxGetN(prhs[0]);
    
    int * NViF, *VCount;
    
    int scalar_dims[2];
    int * FI;
    
    
    #ifdef DEBUG3
    
    printf("NumOfPoints : %i, Num of Faces: %i, Num of Vertices: %i, Num of Max Nbors: %i, Num of Max Vert Faces: %i\n",
            numOfPoints, numOfFaces, numOfVertices, numOfNbors, numOfVertFaces);
    fflush(stdout);
    #endif
    
    if (mxGetM(prhs[0]) != 3)	mexErrMsgTxt("Point should be three dimensional!");
    if (mxGetN(prhs[1]) != mxGetN(prhs[3]) 
        && mxGetN(prhs[1]) != mxGetN(prhs[4])) 
        mexErrMsgTxt("Inconsistent number of vertices!");
    if (mxGetN(prhs[0]) != mxGetN(prhs[5]))
        mexErrMsgTxt("Numof points and number of seedVertices should be the same!");
    
    scalar_dims[0] = 1;
    scalar_dims[1] = numOfPoints;
    plhs[0] = mxCreateNumericArray(2,scalar_dims,mxINT32_CLASS, mxREAL);
    
    FI = (int *) mxGetData(plhs[0]);
    
    
    plhs[1] = mxCreateNumericArray(2, scalar_dims, mxINT32_CLASS, mxREAL);
    
    NViF = (int *) mxGetData(plhs[1]);
       
    MARS_findFaces(points, numOfPoints, vertices, numOfVertices, faces, numOfFaces,
                  vertexNbors, numOfNbors, vertexFaces, numOfVertFaces, seedVertices, FI, NViF);
 
    return;
};
