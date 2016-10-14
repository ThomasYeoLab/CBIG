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

	File:	MARS_findFaces.h
	
    MARS_findFaces( points, numOfPoints, 
                    vertices, numOfVertices,
                    faces, numOfFaces,
                    vertexNbors, numOfMaxNbors, 
                    vertexFaces, numOfMaxFaces,
                    seedVertices,FaceInds_Cstyle )
	
 *   assuming memory has been allocated!!!!
	 point = 3 x numOfPoints (single precision)
     mesh.vertices = 3 x numOfVertices (single) 
 *   mesh.faces    = 3 x numOfFaces    (int32)      (Matlab-style Index)
	 mesh.vertexNbors = N x numOfVertices (int32)   (Matlab-style Index)
 *   mesh.vertexFaces = M x numOfVertices (int32)   (Matlab-style Index)
     seedVertices = 1 x numOfPoints (int32)         (Matlab-style Index)
     FaceInds_CStyle = 1 x numOfPoints (int32)      (Matlab-style index)
      
	 Date:	October 2006

	Copyright (c) 2006 by Mert Rory Sabuncu
************************************************************************/

/* 
 *Local Routines
 */
#include <time.h>
#include <sys/timeb.h>
#include <math.h>
#include "MARS_vec3D.h"

#ifndef MARS_FINDFACES_H
#define MARS_FINDFACES_H


/*#define VERBOSE*/
/*#define MIN_THRESH  -0.000001f*/
#define MIN_THRESH  -0.0001f
#define VERBOSITY 0

#if defined(WIN32) || defined(WIN64)
    typedef __int64 int64;
    typedef unsigned __int64 uint64;
#else
    typedef long long int64;
    typedef unsigned long long uint64;
#endif


static bool isInTriangle(const float * point, 
                            const float * v0, const float * v1, const float * v2, float thresh);

static bool isOnSameSide(const float * p0, const float * p1, const float *a, const float *b, float thresh);

static int  findNextSeed(const float *point, 
                            const float * vertices, const int * vertexNbors, int curSeed ,
                            const bool * isVertexVisited, int numOfMaxNbors, float radius, int numOfVertices);

static int findNearestUnvisitedNbor(const float *point, 
                            const float * vertices, const int * vertexNbors, int curVert ,
                            const bool * isVertexVisited, int numOfMaxNbors, float radius, int numOfVertices);

static bool safetyNet(const float * point, const float * v0, const float * v1, const float * v2);



int64 msTime()
{
#if defined(WIN32) || defined(WIN64)
  struct _timeb tb;
  _ftime(&tb);
#else
  struct timeb tb;
  ftime(&tb);
#endif

  return (int64)((int64)tb.time*(int64)1000) + (int64)tb.millitm;
}


/*
 * The main routine.  
 */
void MARS_findFaces( const float * points, int numOfPoints, 
                const float * vertices, int numOfVertices,
                const int * faces, int numOfFaces,
                const int * vertexNbors, int numOfNbors, 
                const int * vertexFaces, int numOfVertFaces,
                const int * seedVertices, int * FaceInds_MatlabStyle, int * NViF)
{

    /*int numTriTested = 0;*/
    
    bool    *       isVertexVisited, *isFaceVisited;
    int * visitedVertexIndex, * visitedFaceIndex;
    int visitedVertexIndexCnt = 0;
    int visitedFaceIndexCnt = 0;
   
    int *   prev_seeds;
    const float *   cur_point;
    bool            found = false;
    float * v0, *v1, *v2, *v_tmp;
    int f, cur_face_ind, FaceInd;
    float tmp_dist, min_dist, radius, cur_thresh;
    int curSeed, next_seed;
    int p, tmp_cnt,vv;
    int num_fail_safety_net = 0;
    int num_of_trials = 0;
    int max_num_of_trials = 100000;
    int num_of_vertices_w_max_num_of_trials = 0;
    #ifdef DEBUG2
    
    printf("NumOfPoints : %i, Num of Faces: %i, Num of Vertices: %i, Num of Max Nbors: %i, Num of Max Vert Faces: %i\n",
            numOfPoints, numOfFaces, numOfVertices, numOfNbors, numOfVertFaces);
    fflush(stdout);
    #endif
    
    radius = magnitudeVector(vertices);
    
    if (!(isVertexVisited = (bool *) calloc(numOfVertices, sizeof(bool))))
    {
        mexErrMsgTxt("Ran out of memory!");
    }
    if (!(isFaceVisited = (bool *) calloc(numOfFaces, sizeof(bool))))
    {
        mexErrMsgTxt("Ran out of memory!");
    }
    
    if (!(prev_seeds = (int *) calloc(numOfVertices, sizeof(int))))
    {
        mexErrMsgTxt("Ran out of memory!");
    }
    
    if (!(visitedVertexIndex = (int *) calloc(numOfVertices, sizeof(int))))
    {
        mexErrMsgTxt("Ran out of memory!");
    }
    if (!(visitedFaceIndex = (int *) calloc(numOfFaces, sizeof(int))))
    {
        mexErrMsgTxt("Ran out of memory!");
    }
    
    
    /* initialize prev_seeds */
    
    for (vv = 0; vv < numOfVertices; vv++)
        prev_seeds[vv] = -1;
    
    for (p = 0; p < numOfPoints; p++)
    {
        found = false;
        cur_point = &points[3 * p];
        cur_thresh = MIN_THRESH;
        
        #ifdef DEBUG2
        
        if (!(p%100))
            printf("Currently finding face of %i th point\n", p + 1);
        
         fflush(stdout);
        #endif
        
        /*
        * Convert input seed index to C-style index
        */
        
        curSeed = seedVertices[p] - 1;
        
        /*printf("distance is %f, tolerance is %f\n", distanceVector(cur_point, &vertices[3*curSeed]), MIN_TOL);*/
        if (distanceVector(cur_point, &vertices[3*curSeed]) < MIN_TOL)
        {
            min_dist = 0.0;
            /*
             * Note: returning Matlab-style index!
             */
            NViF[p] = curSeed + 1; 
            FaceInd = vertexFaces[numOfVertFaces * curSeed];
            FaceInds_MatlabStyle[p] = FaceInd;  
            found = true;
            continue;
        }
        
        prev_seeds[curSeed] = -1;
        
        tmp_cnt = 0;
        num_of_trials = 0;
        while((!found) && (num_of_trials < max_num_of_trials))
        {   
            num_of_trials++;
            tmp_cnt++;
            f = 0;
            isVertexVisited[curSeed]    = true;
            visitedVertexIndex[visitedVertexIndexCnt] = curSeed;
            visitedVertexIndexCnt++;
            
            #ifdef DEBUG
            printf("Visiting vertex index (Matlab): %i\n", curSeed + 1);
            printf("Vertex coordinates: %f %f %f \n", vertices[3*curSeed], vertices[3*curSeed+ 1], vertices[3*curSeed + 2]);
            #endif
            
            while((!found) 
                    && (f < numOfVertFaces) && (vertexFaces[numOfVertFaces * curSeed + f] != 0))
            {
                /*
                *Convert to C-style index
                */
                
                cur_face_ind = vertexFaces[numOfVertFaces * curSeed + f] - 1;
                #ifdef DEBUG
            
                printf("Checking face # %i\n", cur_face_ind);
                 fflush(stdout);
            
                #endif
            
                if (!isFaceVisited[cur_face_ind])
                {
                    isFaceVisited[cur_face_ind] = true;
                    visitedFaceIndex[visitedFaceIndexCnt] = cur_face_ind;
                    visitedFaceIndexCnt++;
                    
                    /*Note: faces: Matlab-style index - so have to convert to C-style index!*/
                    
            
                    v0  =   (float*) &vertices[3*(faces[3*cur_face_ind] - 1)];
                    v1  =   (float*) &vertices[3*(faces[3*cur_face_ind + 1] - 1)];
                    v2  =   (float*) &vertices[3*(faces[3*cur_face_ind + 2] - 1)];           
                
                    /*
                    *isInTriangle() prefers that v0 is the curSeed - so let's switch!
                    *
                    */
                    if (faces[3*cur_face_ind + 1] == curSeed + 1)
                    {
                        v_tmp = v0;
                        v0  = v1;
                        v1  = v_tmp;
                    }
                    else
                        if(faces[3*cur_face_ind + 2] == curSeed + 1)
                    {
                        v_tmp = v0;
                        v0  = v2;
                        v2  = v_tmp;
                    }
                               
                    found = isInTriangle(cur_point, v0, v1, v2, cur_thresh);
                    /*numTriTested++;*/
                    
                    #ifdef DEBUG
            
                    if (found) printf("Triangle found!\n");
                     fflush(stdout);
            
                    #endif
                };               
                f++;
            }
            
            /*
             * Check whether the found triangle makes sense -- a vertex of that face should be sufficiently close to the point
             * Update threshold
             */
            
            if (found)
            {
                if (distanceVector(cur_point, v0) > 0.25f * radius * sqrt(sqrt(163842.0/numOfVertices)))
                {
                    #ifdef VERBOSE
                    printf("Ooops... The face found was too far away to be the actual face -- possibly on other side of sphere!\n");
                    printf("Changing threshold...\n");
                    #endif
                    found = false;
                    
                    cur_thresh *= 10.0f;
                    
                    /* FREE isVertexVisited (unvisit vertices)*/
                    
                    for (vv = 0; vv < visitedVertexIndexCnt; vv++)
                    {
                        isVertexVisited[visitedVertexIndex[vv]]    = false;
                        prev_seeds[visitedVertexIndex[vv]] = -1;
                        visitedVertexIndex[vv] = 0;
                    }                    
                    visitedVertexIndexCnt = 0;
                    
                                      
                    /* FREE isFaceVisited (unvisit faces)*/
                    
                    for (vv = 0; vv < visitedFaceIndexCnt; vv++)
                    {
                        isFaceVisited[visitedFaceIndex[vv]]    = false;
                        visitedFaceIndex[vv] = 0;
                    }                    
                    visitedFaceIndexCnt = 0;
                   
                    curSeed = seedVertices[p] - 1;
                                       
                    prev_seeds[curSeed] = -1;
                }
                
            }
            /*
            * If found record the face index (Matlab-style)
            */
            if (found)
            {
                /*
                 *
                 */
                v0  =   (float*) &vertices[3*(faces[3*cur_face_ind] - 1)];
                v1  =   (float*) &vertices[3*(faces[3*cur_face_ind + 1] - 1)];
                v2  =   (float*) &vertices[3*(faces[3*cur_face_ind + 2] - 1)];
                
                if (cur_thresh < 5.0f * MIN_THRESH)
                {
                    /*
                     * Threshold had to be changed to find face
                     */
                    
                    if (!safetyNet(cur_point, v0, v1, v2))
                    {
                        num_fail_safety_net++;
                     }
                }
                min_dist = (cur_point[0]-v0[0])*(cur_point[0]-v0[0])
                        +(cur_point[1]-v0[1])*(cur_point[1]-v0[1])
                        +(cur_point[2]-v0[2])*(cur_point[2]-v0[2]);
                
                NViF[p] = faces[3*cur_face_ind];
                
                tmp_dist = (cur_point[0]-v1[0])*(cur_point[0]-v1[0])
                        +(cur_point[1]-v1[1])*(cur_point[1]-v1[1])
                        +(cur_point[2]-v1[2])*(cur_point[2]-v1[2]);
                
                if (tmp_dist < min_dist){
                    min_dist = tmp_dist;
                    NViF[p] = faces[3*cur_face_ind + 1];
                }
                
                tmp_dist = (cur_point[0]-v2[0])*(cur_point[0]-v2[0])
                        +(cur_point[1]-v2[1])*(cur_point[1]-v2[1])
                        +(cur_point[2]-v2[2])*(cur_point[2]-v2[2]);
                
                if (tmp_dist < min_dist)
                {
                    min_dist = tmp_dist;
                    NViF[p] = faces[3*cur_face_ind + 2];
                }
                
                FaceInd =   cur_face_ind + 1;  
                break;
            }           
            else{
                #ifdef DEBUG
                if (vertexFaces[numOfVertFaces * curSeed + f] == 0)
            
                    printf("Ran out of faces for seed #%i!!\n", curSeed);
        
                    printf("Moving onto another seed!\n");  
                     fflush(stdout);
        
                #endif
                next_seed = -1;
                
                while (next_seed == -1)
                { 
                    next_seed = findNextSeed(cur_point, vertices, vertexNbors, curSeed, isVertexVisited, numOfNbors, radius, numOfVertices);
                    if (next_seed == -1)
                        curSeed = prev_seeds[curSeed];
                    if (curSeed < 0)
                    {
                        /*
                        printf("Visited %i vertices\n", tmp_cnt); 
                        mexErrMsgTxt("SCREW YOU!");
                         *
                         */
                        #ifdef VERBOSE
                        printf("Ooops... No face was found with this threshold!\n");
                        printf("Changing threshold...\n");
                        #endif
                        
                        found = false;
                        
                        cur_thresh *= 10.0f;
                        
                        
                        /* FREE isVertexVisited (unvisit vertices)*/
                        
                        for (vv = 0; vv < visitedVertexIndexCnt; vv++)
                        {
                            isVertexVisited[visitedVertexIndex[vv]]    = false;
                            prev_seeds[visitedVertexIndex[vv]] = -1;
                            visitedVertexIndex[vv] = 0;
                        }
                        visitedVertexIndexCnt = 0;
                        
                        
                        /* FREE isFaceVisited (unvisit faces)*/
                        
                        for (vv = 0; vv < visitedFaceIndexCnt; vv++)
                        {
                            isFaceVisited[visitedFaceIndex[vv]]    = false;
                            visitedFaceIndex[vv] = 0;
                        }
                        visitedFaceIndexCnt = 0;
                        
                        next_seed = seedVertices[p] - 1;
                        curSeed = -1;
                        break;
                    }
                    
                }
                prev_seeds[next_seed] = curSeed;
                curSeed = next_seed;
            }
        
            #ifdef DEBUG
        
            printf("Now visiting faces of vertex %i, \n", curSeed);
            fflush(stdout);
            
            #endif
            
        }
        
        if (num_of_trials >= max_num_of_trials)
        {
            num_of_vertices_w_max_num_of_trials++;
            curSeed = seedVertices[p] - 1;
            NViF[p] = curSeed + 1; 
            FaceInd = vertexFaces[numOfVertFaces * curSeed];
            FaceInds_MatlabStyle[p] = FaceInd;       
        }
        
        #ifdef DEBUG
        printf("Recording found face index as %i\n", FaceInd);
         fflush(stdout);
        #endif
        FaceInds_MatlabStyle[p] = FaceInd;
        
        /* FREE isVertexVisited (unvisit vertices)*/
        
        for (vv = 0; vv < visitedVertexIndexCnt; vv++)
        {
            isVertexVisited[visitedVertexIndex[vv]]    = false;
            prev_seeds[visitedVertexIndex[vv]] = -1;
            visitedVertexIndex[vv] = 0;
        }
        visitedVertexIndexCnt = 0;
        
        
        /* FREE isFaceVisited (unvisit faces)*/
        
        for (vv = 0; vv < visitedFaceIndexCnt; vv++)
        {
            isFaceVisited[visitedFaceIndex[vv]]    = false;
            visitedFaceIndex[vv] = 0;
        }
        visitedFaceIndexCnt = 0;

    }; 
    
    if (num_of_vertices_w_max_num_of_trials > 0)
    {
        printf("%i number of vertices exceed the max number of trials\n", num_of_vertices_w_max_num_of_trials);
        fflush(stdout);
    }
    
    if(num_fail_safety_net > 0)
    {
        if(VERBOSITY)
            printf("%d vertices fail safety net!\n", num_fail_safety_net);
    }
    
    if (isVertexVisited)
    {
        free(isVertexVisited);
        isVertexVisited = NULL;
    }
    if (isFaceVisited)
    {
        free(isFaceVisited);
        isFaceVisited = NULL;
    }
    if (prev_seeds)
    {
        free(prev_seeds);
        prev_seeds = NULL;
    }
    if (visitedVertexIndex)
    {
        free(visitedVertexIndex);
        visitedVertexIndex = NULL;
    }
    if (visitedFaceIndex)
    {
        free(visitedFaceIndex);
        visitedFaceIndex = NULL;
    }
      
    /*printf("Tested %d triangles\n", numTriTested);*/
    return;
}

static bool isInTriangle(const float * point, 
                            const float * v0, const float * v1, const float * v2, float thresh)
{
   /*
    *Project point onto plane
    */
    float projected_point[3];
    
    projectPoint(point, v0, v1, v2, projected_point);
          
    return (isOnSameSide(projected_point, v0, v1, v2, thresh) 
            && isOnSameSide(projected_point, v1, v2, v0, thresh)
            && isOnSameSide(projected_point, v2, v0, v1, thresh));
}


    
static bool isOnSameSide(const float * p0, const float * p1, const float *a, const float *b, float thresh)
{
 
    float   edge_ab[3], vec1[3], vec2[3], tmp0[3], tmp1[3];
    
    edge_ab[0] = b[0] - a[0];
    edge_ab[1] = b[1] - a[1];
    edge_ab[2] = b[2] - a[2];
    
    tmp0[0] = p0[0] - a[0];
    tmp0[1] = p0[1] - a[1];
    tmp0[2] = p0[2] - a[2];
    
    tmp1[0] = p1[0] - a[0];
    tmp1[1] = p1[1] - a[1];
    tmp1[2] = p1[2] - a[2];
    
    normalizeVector(tmp0);
    normalizeVector(tmp1);
    normalizeVector(edge_ab);
    
    crossVectors(tmp0, edge_ab, vec1);
    crossVectors(tmp1, edge_ab, vec2);   
    
    return (dotVectors(vec1,vec2) >= thresh);
}

static int  findNextSeed(const float *point, 
                         const float * vertices, const int * vertexNbors, int curSeed ,
                         const bool * isVertexVisited, int numOfMaxNbors, float radius, int numOfVertices)
{

    /*
     *NOTE: input vertexNbors is Matlab-style index!!!
     *AND   curSeed is C-style index!!!
     */
    
    int next_seed = -1;
    next_seed = findNearestUnvisitedNbor(point, vertices,vertexNbors,curSeed,isVertexVisited,numOfMaxNbors, radius, numOfVertices);
    return next_seed;
}

static int findNearestUnvisitedNbor(const float *point, 
                            const float * vertices, const int * vertexNbors, int curVert ,
                            const bool * isVertexVisited, int numOfMaxNbors, float radius, int numOfVertices)
                            
{
    
    /*
     *NOTE: input vertexNbors is Matlab-style index!!!
     *AND   curSeed is C-style index!!!
     */
    
    float min_distance_sq = MAX_DISTANCE_SQ;
    float tmp;
    int i, nbor_ind, nearest_unvisitedNbor;
    i = 0;
    nearest_unvisitedNbor = -1;
    #ifdef DEBUG
        
        printf("Finding next seed from %i (C-style)\n", curVert);
       
    #endif
    while ((i < numOfMaxNbors))
    {
        /*
         *Convert to C-style index
         */
        nbor_ind = vertexNbors[curVert * numOfMaxNbors + i] - 1;
        
        if (nbor_ind >= 0)
        {
           if (!isVertexVisited[nbor_ind])
           {
                tmp = (point[0] - vertices[3*nbor_ind])*(point[0] - vertices[3*nbor_ind])
                    +(point[1] - vertices[3*nbor_ind + 1])*(point[1] - vertices[3*nbor_ind + 1])
                    +(point[2] - vertices[3*nbor_ind + 2])*(point[2] - vertices[3*nbor_ind + 2]);
                if (tmp < min_distance_sq)
                {
                    nearest_unvisitedNbor = nbor_ind;
                    min_distance_sq = tmp;
                }
            }
        }
        else break;
        
        i++;
    }
    
    if (min_distance_sq > (0.25f*radius) * (0.25f*radius) * sqrt(163842.0/numOfVertices))
        nearest_unvisitedNbor = -1;
    
     
    return nearest_unvisitedNbor;
};

static bool safetyNet(const float * point, const float * v0, const float * v1, const float * v2)
{
    float max0, max1, max2, min0, min1, min2;
    float pp[3];
    bool safetyBool;
    
    max0 = min0 = v0[0];
    max1 = min1 = v0[1];
    max2 = min2 = v0[2];
    
    if (v1[0] > max0)
        max0 = v1[0];
    if (v2[0] > max0)
        max0 = v2[0];
    
    if (v1[1] > max1)
        max1 = v1[1];
    if (v2[1] > max1)
        max1 = v2[1];
    
    if (v1[2] > max2)
        max2 = v1[2];
    if (v2[2] > max2)
        max2 = v2[2];
    
    if (v1[0] < min0)
        min0 = v1[0];
    if (v2[0] < min0)
        min0 = v2[0];
    
    if (v1[1] < min1)
        min1 = v1[1];
    if (v2[1] < min1)
        min1 = v2[1];
    
    if (v1[2] < min2)
        min2 = v1[2];
    if (v2[2] < min2)
        min2 = v2[2];
    projectPoint(point, v0, v1, v2, pp);
    safetyBool = pp[0] >= min0 && pp[0] <= max0 && pp[1] >= min1 && pp[1] <= max1 && pp[2] >= min2 && pp[2] <= max2;
    
    if(!safetyBool)
    {       
        if(VERBOSITY)
        {
        printf("v0      : %f, %f, %f\n", v0[0], v0[1], v0[2]);
        printf("v1      : %f, %f, %f\n", v1[0], v1[1], v1[2]);
        printf("v2      : %f, %f, %f\n", v2[0], v2[1], v2[2]);
        printf("point   : %f, %f, %f\n", point[0], point[1], point[2]);
        printf("pp      : %f, %f, %f\n", pp[0], pp[1], pp[2]);
        }
    }
    return (safetyBool);    
    
}
#endif
