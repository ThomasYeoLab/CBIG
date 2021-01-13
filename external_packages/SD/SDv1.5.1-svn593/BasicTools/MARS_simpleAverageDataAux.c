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

	File:   MARS_simpleAverageDataAux.c   
	A mex file for Matlab (c)
	The calling syntax is:
			
	averaged_data = MARS_simpleAverageDataAux(mesh.vertices, mesh.vertexNbors, data, sigma_sq); 
    
	
	 mesh.vertices = 3 x numOfVertices (single) 
 *   mesh.vertexNbors = N x numOfVertices (int32)
 *   data = d x NumOfVertices (single)
 *   sigma_sq = scalar (single)
     averaged_data = d x NumOfVertices (single)
      
	 Date:	October 2006

	Copyright (c) 2006 by Mert Rory Sabuncu
************************************************************************/

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mex.h"
/*
 #define DEBUG
 */

void
mexFunction(
			int nlhs, 
			mxArray *plhs[], 
			int nrhs, 
			const mxArray *prhs[])
{
    const float *   vertices		=	(float*) mxGetPr(prhs[0]);
	const int	*	vertexNbors	    =   (int*)  mxGetPr(prhs[1]);
    const float *   data            =   (float*) mxGetPr(prhs[2]);
    const float     var        =   (float) mxGetScalar(prhs[3]);
    
    const int       numOfVertices   =   mxGetN(prhs[0]);
    const int       numOfVertNbors  =   mxGetM(prhs[1]);
	const int       dataDim         =   mxGetM(prhs[2]);

    int d, nbor_ind, i, n, dd;
    const float * cur_point;
    const float * cur_nbor;
    float * out_data;
    float dist2nbor_sq, w_tot, sc, w_cur;
    float inv_two_var;
    float * tmp_data;
    mwSize out_dims[2];
    
    out_dims[0] = dataDim;
    out_dims[1] = numOfVertices;
    
    if (mxGetN(prhs[0]) != mxGetN(prhs[1]))  mexErrMsgTxt("Vertices and vertexNbors should have the same size[2]!");
    if (mxGetM(prhs[0]) != 3)   mexErrMsgTxt("Number of dimensions should be 3!");
    /*if (mxGetM(prhs[2]) != mxGetM(prhs[3])) mexErrMsgTxt("Data and sigma should be same dimension!");*/
    if (mxGetN(prhs[0]) != mxGetN(prhs[2])) mexErrMsgTxt("data not right size!");
    
    plhs[0] = mxCreateNumericArray(2, out_dims, mxSINGLE_CLASS, mxREAL);
    
    out_data = (float *) mxGetData(plhs[0]);
    
    if (!(tmp_data = (float*) calloc(dataDim, sizeof(float))))
        mexErrMsgTxt("Memory allocation error!");
    
    #ifdef DEBUG
    
    printf("Num of vertices: %i, Num of VertexNbors: %i, Data dimension: %i\n", numOfVertices, numOfVertNbors, dataDim);
    
    #endif
    
    if (var > 0)     
        inv_two_var = 1.0/2.0/var;
    
    
    for (i = 0; i < numOfVertices; i++)
    {
        w_tot = 1;
        for (d = 0; d < dataDim; d++)
            tmp_data[d] = data[dataDim*i + d];
        
        n = 0;
        cur_point =  &vertices[i  * 3];
        
        while (n < numOfVertNbors)
        {
                    /*
                     *Convert to C-style index!
                     */
            nbor_ind = vertexNbors[numOfVertNbors * i + n] - 1;
            
            if (nbor_ind < 0)
                break;
            
            cur_nbor = &vertices[3*nbor_ind];
            dist2nbor_sq = 0.0;
            for (dd = 0; dd < 3; dd++)
                dist2nbor_sq += (cur_point[dd] - cur_nbor[dd]) * (cur_point[dd] - cur_nbor[dd]);
            
            if(var > 0)
                w_cur = (float) (expf(-dist2nbor_sq*inv_two_var));
            else
                w_cur = 1.0;
            
            for (d = 0; d < dataDim; d++)
                tmp_data[d] += w_cur * data[dataDim * nbor_ind + d];
            
            w_tot += w_cur;
            n++;
        }
        for (d = 0; d < dataDim; d++)
        {
            tmp_data[d] = tmp_data[d]/w_tot;
            out_data[dataDim * i + d] = tmp_data[d];
        }
    }

    free(tmp_data);
    return;
}
