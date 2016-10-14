% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================
mex MARS_calculateSurfaceArea.c
mex MARS_convertFaces2FacesOfVert.c
mex MARS_convertFaces2VertNbors.c
mex MARS_findFaceAux.c
mex MARS_findNVAux.c
mex MARS_linearInterpolateAux.c
mex MARS_linearInterpolateAuxWGrad.c
mex MARS_simpleAverageDataAux.c
mex MARS_computeVertexDistSq2Nbors.c
mex MARS_computeFoldingGradAux.c
mex MARS_computeFoldingEnergyAux.c
mex MARS_computeFoldingGradFastAux2.c
mex MARS_computeDiffDataVertex2Nbors.c
mex MARS_computeFoldingGradFastAux.c
mex MARS_computeFoldingEnergyFastAux.c
mex MARS_computeMeshFaceAreas.c
mex MARS_computeAreaEnergyAux.c
mex MARS_computeAreaGradAux.c
mex FindNeighborhoodGivenRadiusAux.cpp
mex MARS_linearInterpolateVertexAuxWGrad.c
mex MARS_linearInterpolateVertexDisplacementWGrad.c

if(ispc)
    check_dir = 'dir';
else
    check_dir = 'ls';
end
[status, y] = system([check_dir fullfile(' ..', 'min_heap')]);
if(status == 0)
    
%     cd('../min_heap');
%     min_heap_compile;
%     cd('../BasicTools');
    
    mex('MARS_DT_Boundary.c', ...
    ['-I' fullfile('..','min_heap')],...
    ['-L' fullfile('..','min_heap')],...
    '-lmin_heap');

    mex('MARS_distSrcs2Dests.c', ...
    ['-I' fullfile('..','min_heap')],...
    ['-L' fullfile('..','min_heap')],...
    '-lmin_heap');

    mex('MARS_src2DestsWithinRange.c', ...
    ['-I' fullfile('..','min_heap')],...
    ['-L' fullfile('..','min_heap')],...
    '-lmin_heap');
else
    disp('min_heap not found, skip compilation of distance transform functions which needs min heap (This is ok)');
end