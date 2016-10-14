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

	File:	MARS_findNVAux.c
	A mex file for Matlab (c)
	The calling syntax is:
			
			 [NVs, distances] = MARS_findNVAux(points, mesh.vertices, mesh.vertexNbors, seedVertices)
	
	 point = 3 x numOfPoints (single precision)
     mesh.vertices = 3 x numOfVertices (single) 
	 mesh.vertexNbors = N x numOfVertices (int32)
     seedVertex = scalar (int32)
     NVs = 1 x numOfPoints (int32)
 *   distances = 1 x numOfPoints (int32)
 
	 Date:	October 2006

	Copyright (c) 2006 by Mert Rory Sabuncu
 *
 *NOTE: This algorithm is NOT guaranteed to find the global nearest vertex in the mesh - especially if it's an irregular mesh!!
 *But should get pretty close!


************************************************************************/

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mex.h"
#include "matrix.h"

#define MAX_DISTANCE_SQ 100000
#define MAX_NUM_RAND 1e3

/*#define DEBUG3
  #define DEBUG*/

/*
 
 *Local Routines
 
 */

static void findNearestVertexAmongstNbors(const float * point, 
                                            const float * vertices, 
                                            const int * vertexNbors, 
                                            int curSeed,
                                            int numOfNbors,
                                            int * nearest, 
                                            float * distance);
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
	
	const float *	points   		=	(float*) mxGetPr(prhs[0]);
	const float *   vertices		=	(float*) mxGetPr(prhs[1]);
	const int	*	vertexNbors	    =   (int*) mxGetPr(prhs[2]);
    const int   *   seedVertices    =   (int*) mxGetPr(prhs[3]);
    
    const int       numOfNbors      =   mxGetM(prhs[2]);
    const int       numOfVertices   =   mxGetN(prhs[1]);
    const int       numOfPoints     =   mxGetN(prhs[0]);
    /*
     bool *          isVertexVisited;
    */
    
    /*Convert input Matlab index to C-style index
     */
    
    float           distance, cur_distance;
    int             curSeed;
    int             scalar_dims[2];
    int             nearest, p;
    int *           NV;
    float *         min_dist;
    const float *   cur_point;
    int rand_count;
    
    /* Added by Thomas to take care of worst case where no seed is found */
    int             bestSeed = -1;
    float           best_distance = MAX_DISTANCE_SQ;
    
    float radius;
    
	if (mxGetM(prhs[0]) != 3)	mexErrMsgTxt("Point should be three dimensional!");
    if (mxGetN(prhs[1]) != mxGetN(prhs[2])) mexErrMsgTxt("Inconsistent number of vertices!");
    if (mxGetN(prhs[0]) != mxGetN(prhs[3])) mexErrMsgTxt("Number of poitns should be the same as # of seedVertices!");
    #ifdef DEBUG3
    
    printf("Number of points: %i, Number of vertices: %i, numOfMaxNbors: %i \n", numOfPoints, numOfVertices, numOfNbors);
    
    printf("Vertex #2: %f, %f, %f \n", vertices[3*2], vertices[3*2+1], vertices[3*2+2]);
    
    fflush(stdout);
    
    #endif
    
    /*
    if (!(isVertexVisited = (bool *) calloc(numOfVertices, sizeof(bool)))
        mexErrMsgTxt("Ran out of memory!");
    */
    radius = vertices[0]*vertices[0] + vertices[1]*vertices[1] + vertices[2]*vertices[2];
    
    scalar_dims[0] = 1;
    scalar_dims[1] = numOfPoints;
    plhs[0] = mxCreateNumericArray(2,scalar_dims,mxINT32_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2,scalar_dims,mxSINGLE_CLASS, mxREAL);
    
    NV = (int *) mxGetData(plhs[0]);
    min_dist = (float*) mxGetData(plhs[1]);
    
    for (p = 0; p < numOfPoints; p++)
    {
        cur_point = &points[3 * p];
        rand_count = 0;
        bestSeed = -1;
        best_distance = MAX_DISTANCE_SQ;
        
        /*
        * Convert input seed index to C-style index
        */
        
        if (seedVertices[p] > 0)
            curSeed = seedVertices[p] - 1;
        else
        {
            if (p == 0)
                curSeed = 0;
            else curSeed = NV[p-1] - 1; /* C-style index!!!*/
        };
    
        cur_distance = (cur_point[0] - vertices[3*curSeed])*(cur_point[0] - vertices[3*curSeed])
                        + (cur_point[1] - vertices[3*curSeed+1])*(cur_point[1] - vertices[3*curSeed+1])
                        + (cur_point[2] - vertices[3*curSeed+2])*(cur_point[2] - vertices[3*curSeed+2]);
        
    
        #ifdef DEBUG
        printf("Sphere radius is: %f\n", sqrt(radius));
        fflush(stdout);
        #endif
    
        srand( time(NULL));
        while (1)
        {       
            findNearestVertexAmongstNbors(cur_point, vertices, vertexNbors, curSeed, numOfNbors, &nearest, &distance);
            #ifdef DEBUG
        
            printf("current seed: %i, current distance: %f, nearest neighbor: %i, distance to NN: %f\n", 
                curSeed, cur_distance, nearest, distance);
            fflush(stdout);
            #endif
            
            
            
            if (distance < cur_distance)
            {
                curSeed = nearest;
                cur_distance = distance;
                
                if(distance < best_distance)
                {
                    best_distance = distance;
                    bestSeed = nearest;
                }
            } 
            else
            {
                /*
                * check if found min distance is acceptable
                */
                if (cur_distance < 0.1f * 0.1f * sqrt(163842.0/numOfVertices) * radius) /*Note that the radius is already squared: so it's more like 100^2*/
                {    
                    if(rand_count >= 100)
                    {
                        printf("RANDOM: vertex %d (matlab index), rand_count: %d\n", p+1, rand_count);
                        fflush(stdout);
                    }
                    break;
                }
                else
                    /*
                    * if the found min distance is not acceptable take a random seed
                    */
                {
                    rand_count++;
                    if(rand_count > numOfVertices)
                    {
                        
                        printf("Problem vertex (matlab index): %d, rand_count: %d, bestSeed: %d, best_distance: %f\n", p+1, rand_count, bestSeed, sqrt(best_distance));
                        fflush(stdout);
                        curSeed = bestSeed;
                        cur_distance = best_distance;
                        break;
                        /*mexErrMsgTxt("MAX_NUM_RAND exceeded");*/
                    }
                    /*printf("RANDOM: Couldn't find nearest vertex with input seed - will try again!\n");
                      fflush(stdout);
                      curSeed = rand_count - 1;*/
                    curSeed = (int) rand()%numOfVertices;
                    cur_distance = (cur_point[0] - vertices[3*curSeed])*(cur_point[0] - vertices[3*curSeed])
                        + (cur_point[1] - vertices[3*curSeed+1])*(cur_point[1] - vertices[3*curSeed+1])
                        + (cur_point[2] - vertices[3*curSeed+2])*(cur_point[2] - vertices[3*curSeed+2]);
                    
                
                }            
            }
        }
    
        
        /* Convert to Matlab-style index*/
        curSeed += 1;
        /*Convert to proper Eucledian distance*/
        cur_distance = sqrt(cur_distance);
    
        #ifdef DEBUG
    
        printf("Final seed: %i, final distance: %f: \n", 
                curSeed, cur_distance);
        fflush(stdout);
    
        #endif 
    
        NV[p] = curSeed;
        min_dist[p] = cur_distance;
    }
    
    return;
}
        

static void findNearestVertexAmongstNbors(const float * point, 
                                          const float * vertices, 
                                          const int * vertexNbors, 
                                          int curSeed,
                                          int numOfNbors,
                                          int * p_nearest, 
                                          float * p_min_dist)
{   
    int i, curNbor;
    float dist;
    (*p_min_dist) = MAX_DISTANCE_SQ;
    
    for (i = 0; i < numOfNbors; i++)
    {
        curNbor = vertexNbors[curSeed * numOfNbors + i] - 1;
        
        if (curNbor >= 0)
        {
            dist = (point[0] - vertices[3*curNbor])*(point[0] - vertices[3*curNbor])
                        + (point[1] - vertices[3*curNbor+1])*(point[1] - vertices[3*curNbor+1])
                        + (point[2] - vertices[3*curNbor+2])*(point[2] - vertices[3*curNbor+2]);
            
            if (dist < (*p_min_dist))
            {
                (*p_min_dist) = dist;
                (*p_nearest) = curNbor;
            }
        }
        
        else 
            break;    
    }
    return;  
};


