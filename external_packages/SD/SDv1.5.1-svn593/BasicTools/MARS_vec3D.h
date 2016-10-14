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
/*
MARS_vec3D.h
*/

#ifndef MARS_VEC3D_H
#define MARS_VEC3D_H

static void crossVectors(const float * u, const float * v, float * res)
{
     /*
     *w = [u(2)*v(3) - u(3)*v(2); u(3)*v(1) - u(1)*v(3); u(1)*v(2) - u(2)*v(1)];
     */
    res[0] = u[1]*v[2] - u[2]*v[1];
    res[1] = u[2]*v[0] - u[0]*v[2];
    res[2] = u[0]*v[1] - u[1]*v[0];
}

static void crossVectors_double(const double * u, const double * v, double * res)
{
     /*
     *w = [u(2)*v(3) - u(3)*v(2); u(3)*v(1) - u(1)*v(3); u(1)*v(2) - u(2)*v(1)];
     */
    res[0] = u[1]*v[2] - u[2]*v[1];
    res[1] = u[2]*v[0] - u[0]*v[2];
    res[2] = u[0]*v[1] - u[1]*v[0];
}

static void crossVectorsNorm(const float * u, const float * v, float * res)
{
    
    float   norm_res; 
    crossVectors(u,v,res);
    
    norm_res = (float) sqrt((float)(res[0] * res[0] + res[1] * res[1] + res[2] * res[2]));
    
    if (norm_res == 0){  
        mexWarnMsgTxt("Normal equal to zero! But should be ok.");
    }
    else 
    {
        res[0] = res[0]/norm_res;
        res[1] = res[1]/norm_res;
        res[2] = res[2]/norm_res;
    }
}

static void crossVectorsNorm_double(const double * u, const double * v, double * res)
{
    
    double   norm_res; 
    crossVectors_double(u,v,res);
    
    norm_res = (double) sqrt((res[0] * res[0] + res[1] * res[1] + res[2] * res[2]));
    
    if (norm_res == 0)  {
        mexWarnMsgTxt("Normal equal to zero! But should be ok.");
    }
    else
    {
        res[0] = res[0]/norm_res;
        res[1] = res[1]/norm_res;
        res[2] = res[2]/norm_res;
    }
}


static float    dotVectors(const float * u, const float * v)
{
    float ans;
    ans = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    
    return ans;
}

static double    dotVectors_double(const double * u, const double * v)
{
    double ans;
    ans = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    
    return ans;
}

static float    magnitudeVector(const float * u)
{
    return (float) (sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

static double    magnitudeVector_double(const double * u)
{
    return (sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]));
}

static void normalizeVector(float * u)
{
 
    float mag = magnitudeVector(u);
    u[0] = u[0]/mag;
    u[1] = u[1]/mag;
    u[2] = u[2]/mag;
    
 }
 
 static void normalizeVector_double(double * u)
{
 
    double mag = magnitudeVector_double(u);
    u[0] = u[0]/mag;
    u[1] = u[1]/mag;
    u[2] = u[2]/mag;
    
 }
 
 static float distanceVector(const float * u, const float *v)
{
    return ((float) (sqrt((u[0] - v[0])*(u[0] - v[0]) + (u[1] - v[1])*(u[1] - v[1]) + (u[2] - v[2])*(u[2] - v[2]))));
}

static double distanceVector_double(const double * u, const double *v)
{
    return ((sqrt((u[0] - v[0])*(u[0] - v[0]) + (u[1] - v[1])*(u[1] - v[1]) + (u[2] - v[2])*(u[2] - v[2]))));
}
/*
static void projectPoint(const float * point, const float * v0, const float * v1, const float * v2, float * projected_point)
{
    float   v0v1[3], v0v2[3], v0point[3], norm[3];
    float   s;
        
    v0v1[0]    =   v1[0] - v0[0];
    v0v1[1]    =   v1[1] - v0[1];
    v0v1[2]    =   v1[2] - v0[2];
    
    v0v2[0]    =   v2[0] - v0[0];
    v0v2[1]    =   v2[1] - v0[1];
    v0v2[2]    =   v2[2] - v0[2];
    
    v0point[0] =   point[0] - v0[0];
    v0point[1] =   point[1] - v0[1];
    v0point[2] =   point[2] - v0[2];
   
    crossVectorsNorm(v0v1, v0v2, norm);
    s = dotVectors(norm, v0point);
    
    projected_point[0] = point[0] - s *norm[0];
    projected_point[1] = point[1] - s *norm[1];
    projected_point[2] = point[2] - s *norm[2];
    return;
}
*/

static void projectPoint(const float * point, const float * v0, const float * v1, const float * v2, float * projected_point)
{
    float   v0v1[3], v0v2[3], norm[3];
    float   s;
        
    v0v1[0]    =   v1[0] - v0[0];
    v0v1[1]    =   v1[1] - v0[1];
    v0v1[2]    =   v1[2] - v0[2];
    
    
    v0v2[0]    =   v2[0] - v0[0];
    v0v2[1]    =   v2[1] - v0[1];
    v0v2[2]    =   v2[2] - v0[2];
    
    
    crossVectorsNorm(v0v1, v0v2, norm);
    if(magnitudeVector(norm) == 0)
    {
        norm[0] = v0[0];
        norm[1] = v0[1];
        norm[2] = v0[2];
        normalizeVector(norm);
    }
    
    s = dotVectors(v0, norm)/dotVectors(point,norm);
    
    projected_point[0] = s*point[0];
    projected_point[1] = s*point[1];
    projected_point[2] = s*point[2];
    return;
}
static void projectPoint_double(const double * point, const double * v0, const double * v1, const double * v2, double * projected_point)
{
    double   v0v1[3], v0v2[3], norm[3];
    double   s;
        
    v0v1[0]    =   v1[0] - v0[0];
    v0v1[1]    =   v1[1] - v0[1];
    v0v1[2]    =   v1[2] - v0[2];
    
    
    v0v2[0]    =   v2[0] - v0[0];
    v0v2[1]    =   v2[1] - v0[1];
    v0v2[2]    =   v2[2] - v0[2];
    
    
    crossVectorsNorm_double(v0v1, v0v2, norm);
    if(magnitudeVector_double(norm) == 0)
    {
        norm[0] = v0[0];
        norm[1] = v0[1];
        norm[2] = v0[2];
        normalizeVector_double(norm);
    }
    
    s = dotVectors_double(v0, norm)/dotVectors_double(point,norm);
    
    projected_point[0] = s*point[0];
    projected_point[1] = s*point[1];
    projected_point[2] = s*point[2];
    return;
}



static void projectPointReturnNorm(const float * point, const float * v0, const float * v1, const float * v2, float * projected_point, float *norm)
{    
    float   v0v1[3], v0v2[3];
    float   s;
        
    v0v1[0]    =   v1[0] - v0[0];
    v0v1[1]    =   v1[1] - v0[1];
    v0v1[2]    =   v1[2] - v0[2];
    
    
    v0v2[0]    =   v2[0] - v0[0];
    v0v2[1]    =   v2[1] - v0[1];
    v0v2[2]    =   v2[2] - v0[2];
    
    
    crossVectorsNorm(v0v1, v0v2, norm);
    if(magnitudeVector(norm) == 0)
    {
        norm[0] = v0[0];
        norm[1] = v0[1];
        norm[2] = v0[2];
        normalizeVector(norm);
    }
    
    s = dotVectors(v0, norm)/dotVectors(point,norm);
    
    projected_point[0] = s*point[0];
    projected_point[1] = s*point[1];
    projected_point[2] = s*point[2];
    return;
}

static void projectPointReturnNorm_double(const double * point, const double * v0, const double * v1, const double * v2, double * projected_point, double *norm)
{    
    double   v0v1[3], v0v2[3];
    double   s;
        
    v0v1[0]    =   v1[0] - v0[0];
    v0v1[1]    =   v1[1] - v0[1];
    v0v1[2]    =   v1[2] - v0[2];
    
    
    v0v2[0]    =   v2[0] - v0[0];
    v0v2[1]    =   v2[1] - v0[1];
    v0v2[2]    =   v2[2] - v0[2];
    
    
    crossVectorsNorm_double(v0v1, v0v2, norm);
    if(magnitudeVector_double(norm) == 0)
    {
        norm[0] = v0[0];
        norm[1] = v0[1];
        norm[2] = v0[2];
        normalizeVector_double(norm);
    }
    
    s = dotVectors_double(v0, norm)/dotVectors_double(point,norm);
    
    projected_point[0] = s*point[0];
    projected_point[1] = s*point[1];
    projected_point[2] = s*point[2];
    return;
}


/*
 * MultiplyMatrix(const float *mat1, const float *mat2, float *out, int dim1, int dim2, int dim3)
 *
 * mat1 is dim1 x dim2
 * mat2 is dim2 x dim3
 * out is dim1 x dim3
 *
 * out = mat1 * mat2
 */

static void MultiplyMatrix(const float *mat1, const float *mat2, float *out, int dim1, int dim2, int dim3)
{
    int i, j, k;
    float total;
    
    for(i = 0; i < dim1; i++){
        for(j = 0; j < dim3; j++){
            
            total = 0;
            for(k = 0; k < dim2; k++)
                total += mat1[i + k*dim1] * mat2[k + j*dim2];
            
            out[i + j*dim1] = total;
            
        }
    }
}




/*
static void projectPointReturnNorm(const float * point, const float * v0, const float * v1, const float * v2, float * projected_point, float *norm)
{
    float   v0v1[3], v0v2[3], v0point[3];
    float   s;
        
    v0v1[0]    =   v1[0] - v0[0];
    v0v1[1]    =   v1[1] - v0[1];
    v0v1[2]    =   v1[2] - v0[2];
    
    v0v2[0]    =   v2[0] - v0[0];
    v0v2[1]    =   v2[1] - v0[1];
    v0v2[2]    =   v2[2] - v0[2];
    
    v0point[0] =   point[0] - v0[0];
    v0point[1] =   point[1] - v0[1];
    v0point[2] =   point[2] - v0[2];
   
    crossVectorsNorm(v0v1, v0v2, norm);
    s = dotVectors(norm, v0point);
    
    projected_point[0] = point[0] - s *norm[0];
    projected_point[1] = point[1] - s *norm[1];
    projected_point[2] = point[2] - s *norm[2];
    return;
}
 */
#endif
