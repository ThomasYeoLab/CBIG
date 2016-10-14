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

	File:	MARS_computeFoldingGradFastAux.c
	A mex file for Matlab (c)
	The calling syntax is:
			
	[energy, grad, vertexList_out] = MARS_computeFoldingGradFastAux(vertices, mesh.faces, mesh.vertexFaces, vertexList_in, AREA_0)
	
	 mesh.vertices = 3 x numOfVertices (single) 
 *   mesh.faces    = 3 x numOfFaces    (int32)
	 mesh.vertexFaces = M x numOfVertices (int32)
 *   AREA_0 = scalar (single)
     
    [energy] = scalar (single)
 *  [grad] = 3 x numOfVertices (single)
 *
 *  NOTE: This function assumes a certain orientation of faces- which is in the regular icosahedral mesh.
 *          
      
	 Date:	November 2006

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
#include "MARS_linearInterp.h"
#define MIN_TOL     0.0000001f

/*
         * This assumes that in the original mesh (that contains no folds) the faces were saved such that they're oriented inwards.
         * This is the cased on the Icosahedral Mesh !!! Other subject meshes seem to be oriented the other way round!
         *
*/
#define INWARD_NORMAL

static bool isFlipped(float area)
{
    
    #ifdef INWARD_NORMAL
    if (area >  MIN_TOL) return true;
    else return false;
    #endif
    
    #ifdef OUTWARD_NORMAL
    if (area <  - MIN_TOL) return true;
    else return false;
    #endif
}



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

    const float *   vertices		=	(float*) mxGetPr(prhs[0]);
    const int   *   faces           =   (int*)  mxGetPr(prhs[1]);
    const int   *   vertexFaces     =   (int*)  mxGetPr(prhs[2]);
    const int   *   vertexList_in   =   (int*)  mxGetPr(prhs[3]);     
    
    const float     area_0          =   (float) mxGetScalar(prhs[4]);
    
    const int       numOfVertFaces  =   mxGetM(prhs[2]);
    const int       numOfVertices   =   mxGetN(prhs[0]);
    const int       numOfFaces      =   mxGetN(prhs[1]);
    const int       size_list_in    =   mxGetM(prhs[3]);

    int f, scalar_dim, v, cur_face, size_list_out, N_in, vert;
    int *   vertexList_out;
    const float * v0, * v1, * v2;
    int v0_ind, v1_ind, v2_ind, tmp_ind;
    float tmp_vec[3], v0v1[3], v0v2[3], tmp_grad[3];
    float * faceArea,  * energy, *gradient;
    int out_grad_dims[2];
    bool *  isVertexListed;
    int*    flippedVertexList;
    out_grad_dims[0] = 3;
    out_grad_dims[1] = numOfVertices;
    
    scalar_dim = 1;
    size_list_out = 0;
    
    if (mxGetM(prhs[0]) != 3) mexErrMsgTxt("Vertices should be 3 dimensional");
    if (mxGetM(prhs[1]) != 3) mexErrMsgTxt("Faces should have three vertices");
    if (mxGetN(prhs[2]) != mxGetN(prhs[0])) mexErrMsgTxt("Vertices and VertexFaces should be the same number!");
    
    if (!(faceArea = (float*) calloc(numOfFaces, sizeof(float))))
        mexErrMsgTxt("Memory allocation error!");
    
    if (!(isVertexListed = (bool*) calloc(numOfVertices, sizeof(bool))))
        mexErrMsgTxt("Memory allocation error!");
    
    if (!(flippedVertexList = (int*) calloc(numOfVertices, sizeof(int))))
        mexErrMsgTxt("Memory allocation error!");
    
    plhs[0] = mxCreateNumericArray(1,&scalar_dim,mxSINGLE_CLASS, mxREAL);
    
    energy = (float *) mxGetData(plhs[0]);
    
    plhs[1] = mxCreateNumericArray(2, out_grad_dims, mxSINGLE_CLASS, mxREAL);
    
    gradient = (float *) mxGetData(plhs[1]);
    
    /*
     * First pass thru faces :  determine folded faces
     *
     */
    if (size_list_in == 0)  
        N_in = numOfVertices;
    else
        N_in = size_list_in;
    for (vert = 0; vert < N_in; vert ++)
    {
        if (size_list_in == 0)
            v = vert;
        else
            /* Convert to C-style index*/
            v = vertexList_in[vert] - 1; 
        
        for (f = 0; f < numOfVertFaces; f++)
        {
            /*
             *Convert to C-style index!
             */
            cur_face = vertexFaces[numOfVertFaces * v + f] - 1;
            
            if (cur_face < 0)
                f = numOfVertFaces;
            else
            {
                if (faceArea[cur_face] == 0.0)
                {   
                    /*
                     *Notice the following are MATLAB-style indices
                     */
                    v0_ind = faces[cur_face * 3];
                    v1_ind = faces[cur_face * 3 + 1];
                    v2_ind = faces[cur_face * 3 + 2];
                    v0 = &vertices[3*(v0_ind - 1)];
                    v1 = &vertices[3*(v1_ind - 1)];
                    v2 = &vertices[3*(v2_ind - 1)];
                    
                    faceArea[cur_face] = orientedArea(v0, v1,v2, v0);
                    if (isFlipped(faceArea[cur_face]))
                    {
                        energy[0] += faceArea[cur_face] + area_0;
                        
                        if (!isVertexListed[v0_ind - 1])
                        {
                            flippedVertexList[size_list_out] = v0_ind;
                            isVertexListed[v0_ind - 1] = true;
                            size_list_out++;
                        }
                        
                        if (!isVertexListed[v1_ind - 1])
                        {
                            flippedVertexList[size_list_out] = v1_ind;
                            isVertexListed[v1_ind - 1] = true;
                            size_list_out++;
                        }
                        
                        if (!isVertexListed[v2_ind - 1])
                        {
                            flippedVertexList[size_list_out] = v2_ind;
                            isVertexListed[v2_ind - 1] = true;
                            size_list_out++;
                        }                        
                        
                    }
                }
            }
        }
    }
    
    plhs[2] = mxCreateNumericArray(1,&size_list_out,mxINT32_CLASS, mxREAL);
    vertexList_out = (int*) mxGetData(plhs[2]);
    
    
    for (vert = 0; vert < size_list_out; vert++)
    {
        
        if (flippedVertexList[vert] <= 0)
        {
            break;
        }
        else
        {
            /*
             *Convert to C-style index!!!
             */
            v = flippedVertexList[vert] - 1;
            /*
             *return Matlab-style index!!
             */
            vertexList_out[vert] = flippedVertexList[vert];
            
            
            for (f = 0; f < numOfVertFaces; f++)
            {
            /*
             *Convert to C-style index!
             */
                cur_face = vertexFaces[numOfVertFaces * v + f] - 1;
                
                if (cur_face < 0)
                    f = numOfVertFaces;
                else
                {
                    if (isFlipped(faceArea[cur_face]))
                    {
                        v0_ind = faces[cur_face * 3] - 1;
                        v1_ind = faces[cur_face * 3 + 1] - 1;
                        v2_ind = faces[cur_face * 3 + 2] - 1;
                        
                    /*
                     * v0 should be v -- for gradient
                     *
                     */
                        
                        if (v0_ind != v)
                        {
                            tmp_ind = v0_ind;
                            v0_ind = v;
                            
                            if (v1_ind == v)
                            {
                                v1_ind = tmp_ind;
                            }
                            if(v2_ind == v)
                            {
                                v2_ind = tmp_ind;
                            }
                        }
                        
                        v0 = &vertices[3*v0_ind];
                        v1 = &vertices[3*v1_ind];
                        v2 = &vertices[3*v2_ind];
                        
                        computeAreaWGradient(v0, v1, v2, tmp_grad);
                        
                        gradient[3 * v] += tmp_grad[0];
                        gradient[3 * v + 1] += tmp_grad[1];
                        gradient[3 * v + 2] += tmp_grad[2];
                    }
                    
                }
                
            }
        }
            
    }
    
    free(faceArea);
    free(flippedVertexList);
    free(isVertexListed);
    return;
}

