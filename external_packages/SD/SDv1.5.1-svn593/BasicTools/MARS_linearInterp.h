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
#ifndef MARS_LINEARINTERP_H
#define MARS_LINEARINTERP_H

/*
 *linearInterp.h
 */
#define MIN_AREA_TOL 0.9999999
#define VAL_PRECISION 0.00001
#include "MARS_vec3D.h"

static float computeArea(const float v0[3], const float v1[3], const float v2[3])
{
    float res[3], v0v1[3], v0v2[3];
    
    v0v1[0] = v1[0] - v0[0];
    v0v1[1] = v1[1] - v0[1];
    v0v1[2] = v1[2] - v0[2];
    
    v0v2[0] = v2[0] - v0[0];
    v0v2[1] = v2[1] - v0[1];
    v0v2[2] = v2[2] - v0[2];
    
    crossVectors(v0v1, v0v2, res);
    
    return ((float)(sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2])/2.0f));
}

static float computeAreaWGradient(const float v0[3], const float v1[3], const float v2[3], float grad[3])
{
    float res[3], area;
    float v1v2[3], v1v0[3], tmp[3];
    float sc;
    
    v1v2[0] = v2[0] - v1[0];
    v1v2[1] = v2[1] - v1[1];
    v1v2[2] = v2[2] - v1[2];
    
    v1v0[0] = v0[0] - v1[0];
    v1v0[1] = v0[1] - v1[1];
    v1v0[2] = v0[2] - v1[2];
      
    crossVectors(v1v2, v1v0, res);
    
    area = (float)(sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2])/2.0f);
    
    crossVectors(v1v2,v1v0,tmp);
    
    #ifdef DEBUG
    printf("Computing area:\n");
    printf("v0: %f, %f, %f\n", v0[0], v0[1], v0[2]);
    printf("v1: %f, %f, %f\n", v1[0], v1[1], v1[2]);
    printf("v2: %f, %f, %f\n", v2[0], v2[1], v2[2]);
    #endif
    
    if (area >= 1.0-MIN_AREA_TOL)
        sc = (float) (0.25f/area);
    else sc = 0.0f;
    
    grad[0] = (float) sc * (tmp[1] * (v2[2] - v1[2]) + tmp[2] * (v1[1] - v2[1]));
    grad[1] = (float) sc * (tmp[0] * (v1[2] - v2[2]) + tmp[2] * (v2[0] - v1[0]));
    grad[2] = (float) sc * (tmp[0] * (v2[1] - v1[1]) + tmp[1] * (v1[0] - v2[0]));  
   
    return area;
    
}
/*
 *orientedArea: returns positive area if the triangle is oriented in the same direction as the normal (n)
 *
 */

static float orientedArea(const float v0[3], const float v1[3], const float v2[3], const float n[3])
{
    float v0v1[3], v0v2[3], tmp[3], area;
    
    v0v1[0] = v1[0] - v0[0];
    v0v1[1] = v1[1] - v0[1];
    v0v1[2] = v1[2] - v0[2];
    
    v0v2[0] = v2[0] - v0[0];
    v0v2[1] = v2[1] - v0[1];
    v0v2[2] = v2[2] - v0[2];
    
    crossVectors(v0v1, v0v2, tmp);
    area = magnitudeVector(tmp)/2.0f;
    
    if (dotVectors(tmp, n) > 0) 
        return area;
    else return -1.0f * area;   
    
}

static void linearInterp(const float point[3], const float v0[3], const float v1[3], const float v2[3], 
                            const float d0[3], const float d1[3], const float d2[3], int data_dim, float new_data[3])
{
    /*
     *Assumes new_data has been allocated - otherwise will crash!!!
     *
     */
    int i = 0;
    int j = 0;
    float A0, A1, A2, projected_point[3],totalA;
    
    if ( new_data == NULL ) 
        mexErrMsgTxt("Memory not allocated for interpolated data!!"); 
    
        
    /*
     *project point onto face-plane
     */
    
    projectPoint(point, v0, v1, v2, projected_point);
    
    A0 = computeArea(projected_point, v1, v2);
    A1 = computeArea(projected_point, v0, v2);
    A2 = computeArea(projected_point, v0, v1);
    
    totalA = A0 + A1 + A2;
    #ifdef DEBUG
    
    printf("Interpolating Areas: %f %f %f\n", A0, A1, A2);
    
    #endif
    
    A0 = A0/totalA;
    A1 = A1/totalA;
    A2 = A2/totalA;
    
    
   
    if (A0 > MIN_AREA_TOL){
        
        A0 = 1.0;
        A1 = 0.0;
        A2 = 0.0;
    }
    if(A1 > MIN_AREA_TOL)
    {
        A0 = 0.0; 
        A1 = 1.0;
        A2 = 0.0;
    }
    if(A2 > MIN_AREA_TOL)
    {
        A0 = 0.0; 
        A1 = 0.0;
        A2 = 1.0;
    }
    
    for (i = 0; i < data_dim; i++)
    {
        new_data[i] = A0 * d0[i] + A1 * d1[i] + A2 * d2[i];
    }
    
};

/*
 * Compute Chain Rule Due to Projection onto Triangle:
 * 
 * <v0, n>/<p, n> I - <v0, n>/<p, n>^2 pn^T = [vec_x vec_y vec_z], column vectors
 */
static void computeChainRuleFromProjOntoTri(const float *v0, const float *point, const float *norm, float *vec_x, float *vec_y, float *vec_z)
{
    float C1, C2;
    
    C1 = dotVectors(v0, norm)/dotVectors(point, norm);
    C2 = C1/dotVectors(point, norm);
    
     /*First handle Rx*/
    vec_x[0] = C1  - C2*point[0]*norm[0];
    vec_x[1] = 0.0 - C2*point[1]*norm[0];
    vec_x[2] = 0.0 - C2*point[2]*norm[0];
    
    /*Now handle Ry*/
    vec_y[0] = 0.0 - C2*point[0]*norm[1];
    vec_y[1] = C1  - C2*point[1]*norm[1];
    vec_y[2] = 0.0 - C2*point[2]*norm[1];
    
    /*Finally Rz*/
    vec_z[0] = 0.0 - C2*point[0]*norm[2];
    vec_z[1] = 0.0 - C2*point[1]*norm[2];
    vec_z[2] = C1  - C2*point[2]*norm[2];
    
}

/*
 * Compute Chain Rule Due to Projection onto Triangle:
 * 
 * Same as computeChainRuleFromProjOntoTri, but output as matrix.
 *
 * <v0, n>/<p, n> I - <v0, n>/<p, n>^2 pn^T = [vec_x vec_y vec_z], column vectors = mat
 */
static void computeChainRuleFromProjOntoTriMat(const float *v0, const float *point, const float *norm, float *mat)
{
    float C1, C2;
    
    C1 = dotVectors(v0, norm)/dotVectors(point, norm);
    C2 = C1/dotVectors(point, norm);
    
     /*First handle Rx*/
    mat[0] = C1  - C2*point[0]*norm[0];
    mat[1] = 0.0 - C2*point[1]*norm[0];
    mat[2] = 0.0 - C2*point[2]*norm[0];
    
    /*Now handle Ry*/
    mat[3] = 0.0 - C2*point[0]*norm[1];
    mat[4] = C1  - C2*point[1]*norm[1];
    mat[5] = 0.0 - C2*point[2]*norm[1];
    
    /*Finally Rz*/
    mat[6] = 0.0 - C2*point[0]*norm[2];
    mat[7] = 0.0 - C2*point[1]*norm[2];
    mat[8] = C1  - C2*point[2]*norm[2];
    
}

/*
 * Compute Chain Rule Due to Projection onto Triangle:
 * 
 * Same as computeChainRuleFromProjOntoTriMat, but assumes v0 = point.
 *
 * <v0, n>/<p, n> I - <v0, n>/<p, n>^2 pn^T = [vec_x vec_y vec_z], column vectors = mat
 *
 */
static void computeChainRuleFromProjOntoTriMatVertex(const float *point, const float *norm, float *mat)
{
    float C;
    
    C = 1.0/dotVectors(point, norm);
    
     /*First handle Rx*/
    mat[0] = 1   - C*point[0]*norm[0];
    mat[1] = 0.0 - C*point[1]*norm[0];
    mat[2] = 0.0 - C*point[2]*norm[0];
    
    /*Now handle Ry*/
    mat[3] = 0.0 - C*point[0]*norm[1];
    mat[4] = 1   - C*point[1]*norm[1];
    mat[5] = 0.0 - C*point[2]*norm[1];
    
    /*Finally Rz*/
    mat[6] = 0.0 - C*point[0]*norm[2];
    mat[7] = 0.0 - C*point[1]*norm[2];
    mat[8] = 1   - C*point[2]*norm[2];
    
}


static void linearInterpWGrad(const float point[3], const float v0[3], const float v1[3], const float v2[3], 
                            const float d0[3], const float d1[3], const float d2[3], int data_dim, float * new_data, float * grad)
{
    /*Assumes new_data has been allocated - otherwise will crash!!!*/
    
    int i = 0;
    int j = 0;
    float A0, A1, A2, projected_point[3], norm[3], totalA, dA0[3], dA1[3], dA2[3];
    float *temp_grad, temp_vec_x[3], temp_vec_y[3], temp_vec_z[3];
    
    if ( new_data == NULL ) 
        mexErrMsgTxt("Memory not allocated for interpolated data!!"); 
    
        
    
    /*project point onto face-plane*/
    projectPointReturnNorm(point, v0, v1, v2, projected_point, norm);
      
    A0 = computeAreaWGradient(projected_point, v1, v2, dA0);
    A1 = computeAreaWGradient(projected_point, v0, v2, dA1);
    A2 = computeAreaWGradient(projected_point, v0, v1, dA2);
    
    
    
    
    /*To be sure they adds up to 0.*/
    dA2[0] = -dA0[0] - dA1[0];
    dA2[1] = -dA0[1] - dA1[1];
    dA2[2] = -dA0[2] - dA1[2];
    
    totalA = A0 + A1 + A2;
    
    A0 = A0/totalA;
    A1 = A1/totalA;
    A2 = A2/totalA;   
    
   
    if (A0 > MIN_AREA_TOL){
        
        A0 = 1.0;
        A1 = 0.0;
        A2 = 0.0;
    }
    if(A1 > MIN_AREA_TOL)
    {
        A0 = 0.0; 
        A1 = 1.0;
        A2 = 0.0;
    }
    if(A2 > MIN_AREA_TOL)
    {
        A0 = 0.0; 
        A1 = 0.0;
        A2 = 1.0;
    }
    
    temp_grad = (float *) calloc(data_dim*3, sizeof(float));
    for (i = 0; i < data_dim; i++)

    {
        #ifdef DEBUG
        
        if (i == 22)
        {
            printf("Label value #23: %f %f %f\n", d0[i], d1[i], d2[i]);
        }
        #endif
        new_data[i] = A0 * d0[i] + A1 * d1[i] + A2 * d2[i];
        
        if (fabs(new_data[i] - d0[i]) < VAL_PRECISION || fabs(new_data[i] - d1[i]) < VAL_PRECISION || fabs(new_data[i] - d2[i]) < VAL_PRECISION)
        {
            for (j = 0; j < 3; j++) /*interpolated point is one of the vertex of the triangles*/
                temp_grad[i + j * data_dim] = 0.0;
        }
        else
        {
            /*currently, temp_grad is d x 3, and corresponds to the rate of change of data wrt projected point.*/
            for (j = 0; j < 3; j++)
                temp_grad[i + j * data_dim] = (dA0[j] * d0[i] + dA1[j] * d1[i] + dA2[j] * d2[i])/totalA;
        }
    }
    
    /*Now we will complete the grad matrix*/
    computeChainRuleFromProjOntoTri(v0, point, norm, temp_vec_x, temp_vec_y, temp_vec_z);
    /*temp_vec_x[0] = 1; temp_vec_x[1] = 0; temp_vec_x[2] = 0;*/
    /*temp_vec_y[0] = 0; temp_vec_y[1] = 1; temp_vec_y[2] = 0; */
    /*temp_vec_z[0] = 0; temp_vec_z[1] = 0; temp_vec_z[2] = 1; */
    
    /* First handle Rx*/
    for (i = 0; i < data_dim; i++)
        grad[i + 0*data_dim] = temp_vec_x[0]*temp_grad[i +  0 * data_dim] + temp_vec_x[1]*temp_grad[i + 1 * data_dim] + temp_vec_x[2]*temp_grad[i + 2 * data_dim];
    
    /* First handle Ry*/
    for (i = 0; i < data_dim; i++)
        grad[i + 1*data_dim] = temp_vec_y[0]*temp_grad[i +  0 * data_dim] + temp_vec_y[1]*temp_grad[i + 1 * data_dim] + temp_vec_y[2]*temp_grad[i + 2 * data_dim];
    
    /* First handle Rz*/
    for (i = 0; i < data_dim; i++)
        grad[i + 2*data_dim] = temp_vec_z[0]*temp_grad[i +  0 * data_dim] + temp_vec_z[1]*temp_grad[i + 1 * data_dim] + temp_vec_z[2]*temp_grad[i + 2 * data_dim];
    
    free(temp_grad);
    
}


/*
 * Find the unit normal to the edge v0v1 of triangel v0,v1,v2. 
 * Assume unit Normal is along plane of triangle and it's in the direction pointing towards the remaining vertex
 */
static void computeNormal2EdgeOfTriangle(const float *v0, const float *v1, const float *v2, float *norm2edge)
{
    float v0v1[3], v0v2[3], norm2triangle[3];
    
    v0v1[0] = v1[0] - v0[0];
    v0v1[1] = v1[1] - v0[1];
    v0v1[2] = v1[2] - v0[2];
    normalizeVector(v0v1);
    
    v0v2[0] = v2[0] - v0[0];
    v0v2[1] = v2[1] - v0[1];
    v0v2[2] = v2[2] - v0[2];
    normalizeVector(v0v2);
    
    crossVectors(v0v1, v0v2, norm2triangle); /*First find normal*/
    normalizeVector(norm2triangle);
    
    crossVectors(v0v1, norm2triangle, norm2edge); /*Then find edge to triangle*/
    normalizeVector(norm2edge);
    
    if(dotVectors(v0v2, norm2edge) < 0)  /*first edge is in wrong direction, flip it.*/
    {
        norm2edge[0] = -norm2edge[0];
        norm2edge[1] = -norm2edge[1];
        norm2edge[2] = -norm2edge[2];
    }    
}

/*
 * Find the gradient of barycentric triangle associated with edge v0v1 (and vertex v2)
 */
static void computeGradientOfBarycentricTriangle(const float *v0, const float *v1, const float *v2, float *grad)
{
    float v0v1[3], norm2edge[3];
    float base;
    
    computeNormal2EdgeOfTriangle(v0, v1, v2, norm2edge);
    
    v0v1[0] = v1[0] - v0[0];
    v0v1[1] = v1[1] - v0[1];
    v0v1[2] = v1[2] - v0[2];
    base = magnitudeVector(v0v1);
    
    grad[0] = 0.5*base*norm2edge[0];
    grad[1] = 0.5*base*norm2edge[1];
    grad[2] = 0.5*base*norm2edge[2];    
} 

/*
 * Compute gradient (including chain rule due to warp and projection) of triangel v0, v1, v2. 
 *
 * grad is data_dim x 3
 *
 * It's the same function as linearInterpWGrad, except it's guaranteed to work if point passed in is one of the vertices.
 */
static void barycentricInterpWGrad(const float point[3], const float v0[3], const float v1[3], const float v2[3],  
                            const float d0[3], const float d1[3], const float d2[3], int data_dim, float * new_data, float *grad)
{
    /*Assumes new_data has been allocated - otherwise will crash!!!*/
    
    int i = 0;
    int j = 0;
    float A0, A1, A2, projected_point[3], norm[3], totalA, dA0[3], dA1[3], dA2[3];
    float *temp_grad, temp_mat[9]; 
    
    if ( new_data == NULL ) 
        mexErrMsgTxt("Memory not allocated for interpolated data!!"); 
       
    /*project point onto face-plane*/
    projectPointReturnNorm(point, v0, v1, v2, projected_point, norm);
    
    
    A0 = computeArea(projected_point, v1, v2);
    A1 = computeArea(projected_point, v0, v2);
    A2 = computeArea(projected_point, v0, v1);    
     
    
    computeGradientOfBarycentricTriangle(v0, v1, v2, dA2);
    computeGradientOfBarycentricTriangle(v2, v0, v1, dA1);
    computeGradientOfBarycentricTriangle(v1, v2, v0, dA0);

    totalA = A0 + A1 + A2;
    
    A0 = A0/totalA;
    A1 = A1/totalA;
    A2 = A2/totalA;   
    
    
    temp_grad = (float *) calloc(data_dim*3, sizeof(float));
    for (i = 0; i < data_dim; i++)

    {
        #ifdef DEBUG
        
        if (i == 22)
        {
            printf("Label value #23: %f %f %f\n", d0[i], d1[i], d2[i]);
        }
        #endif
        new_data[i] = A0 * d0[i] + A1 * d1[i] + A2 * d2[i];
        
        for (j = 0; j < 3; j++)
            temp_grad[i + j * data_dim] = (dA0[j] * d0[i] + dA1[j] * d1[i] + dA2[j] * d2[i])/totalA;

    }
    
    computeChainRuleFromProjOntoTriMatVertex(point, norm, temp_mat);
    MultiplyMatrix(temp_grad, temp_mat, grad, data_dim, 3, 3);
    
    free(temp_grad);
}


/*
 * Compute gradient (including chain rule due to warp and projection) of triangel v0, v1, v2. 
 *
 * grad is data_dim x 3
 *
 * It's the same function as barycentricInterpWGrad, except it does not normalize the area.
 * Note that we assume data_dim = 3!!
 * 
 * Also, we have to have a multiplier in front.
 * 
 */
static void barycentricInterpWGradWONormalization(const float point[3], const float v0[3], const float v1[3], const float v2[3],  
                            const float d0[3], const float d1[3], const float d2[3], int data_dim, float * new_data, float *grad)
{
    /*Assumes new_data has been allocated - otherwise will crash!!!*/
    
    int i = 0;
    int j = 0;
    float A0, A1, A2, projected_point[3], norm[3], totalA, dA0[3], dA1[3], dA2[3];
    float *temp_grad, temp_mat[9], temp_grad2[9];
    
    if ( new_data == NULL ) 
        mexErrMsgTxt("Memory not allocated for interpolated data!!"); 
       
    if(data_dim != 3)
        mexErrMsgTxt("Assume data dim is 3!!"); 
    
    
    /*project point onto face-plane*/
    projectPointReturnNorm(point, v0, v1, v2, projected_point, norm);
    
    
    A0 = computeArea(projected_point, v1, v2);
    A1 = computeArea(projected_point, v0, v2);
    A2 = computeArea(projected_point, v0, v1);    
     
    
    computeGradientOfBarycentricTriangle(v0, v1, v2, dA2);
    computeGradientOfBarycentricTriangle(v2, v0, v1, dA1);
    computeGradientOfBarycentricTriangle(v1, v2, v0, dA0);

    totalA = A0 + A1 + A2;
    
    A0 = A0/totalA;
    A1 = A1/totalA;
    A2 = A2/totalA;   
    
    temp_grad = (float *) calloc(data_dim*3, sizeof(float));
    for (i = 0; i < data_dim; i++)

    {
        #ifdef DEBUG
        
        if (i == 22)
        {
            printf("Label value #23: %f %f %f\n", d0[i], d1[i], d2[i]);
        }
        #endif
        new_data[i] = A0 * d0[i] + A1 * d1[i] + A2 * d2[i];
        
        for (j = 0; j < 3; j++)
            temp_grad[i + j * data_dim] = (dA0[j] * d0[i] + dA1[j] * d1[i] + dA2[j] * d2[i]); /* This is the part that is different from barycentricInterpWGrad because of no division by A.*/

    }
    
    /*Now we will complete the grad matrix*/
    computeChainRuleFromProjOntoTriMatVertex(point, norm, temp_mat);
    MultiplyMatrix(temp_grad, temp_mat, temp_grad2, data_dim, 3, 3);
    
    /* This part is different from barycentricInterpWGrad as well.*/
    MultiplyMatrix(temp_mat, temp_grad2, grad, data_dim, 3, 3);
    
    free(temp_grad);
    
/*     if(point[0] == 0 && point[1] == 0 && point[2] == 100)
     {
         printf("d0: %f %f %f\n", d0[0], d0[1], d0[2]);
         printf("d1: %f %f %f\n", d1[0], d1[1], d1[2]);
         printf("d2: %f %f %f\n", d2[0], d2[1], d2[2]);
         printf("dA0: %f %f %f\n", dA0[0], dA0[1], dA0[2]);
         printf("dA1: %f %f %f\n", dA1[0], dA1[1], dA1[2]);
         printf("dA2: %f %f %f\n", dA2[0], dA2[1], dA2[2]);
         printf("**********\n");
     }*/
}


#endif












