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
#define MAX_DISTANCE_SQ 100000
#define MIN_TOL     0.00005

#include "mex.h"
#include "MARS_findFaces.h"
#include "math.h"
#include "time.h"

static const size_t size = 1<<16;
static const int iteration = 1000000;


void
mexFunction(
			int nlhs, 
			mxArray *plhs[], 
			int nrhs, 
			const mxArray *prhs[])
{
    
    srand( time(NULL));
    
    float point[3];
    float v0[3];
    float v1[3];
    float v2[3];
    int i;
    
    
    
    v0[0] = 20.8820;
    v0[0] = -89.1166;
    v0[1] = -40.2763;
    
    v1[0] = 20.4478;
    v1[0] = -89.3205;
    v1[1] = -40.0466;
    
    v2[0] = 20.9045;
    v2[0] = -88.8135;
    v2[1] = -40.9288;


    for (i = 0; i < iteration; i++)
    {
        point[0] = ((double) rand())/RAND_MAX * 100;
        point[1] = ((double) rand())/RAND_MAX * 100;
        point[2] = ((double) rand())/RAND_MAX * 100;
        
        isInTriangle(point, v0,  v1, v2, 1e-5);
    
    }
    
}
