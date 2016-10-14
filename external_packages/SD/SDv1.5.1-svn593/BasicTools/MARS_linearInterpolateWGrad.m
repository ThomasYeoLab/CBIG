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
function [vals, grads, NViF, NF] = MARS_linearInterpolateWGrad(points, mesh, data, nearest_vertex_guess)
%function [vals] = MARS_linearInterpolate(point, mesh, data, nearest_vertex_guess)

%description: Returns interpolated value(s) of data (on mesh) at point

%MARS Project (Thomas Yeo and Mert Sabuncu), MIT, CSAIL, (c) 2006
%Author: Mert

if (nargin < 4)
    [nearest_vertex_guess, d] =  MARS_findNV_kdTree(single(points), mesh.vertices);
    %nearest_vertex_guess =  MARS_findNearestVertex(single(points), mesh);
else 
    [nearest_vertex_guess, d] =  MARS_findNV_kdTree(single(points), mesh.vertices);
    %nearest_vertex_guess =  MARS_findNearestVertex(single(points), mesh, nearest_vertex_guess);
end

[vals, grads, NViF, NF] = MARS_linearInterpolateAuxWGrad(single(points), mesh.vertices, mesh.faces, mesh.vertexNbors,... 
                                      mesh.vertexFaces, nearest_vertex_guess, single(data));

if(~isempty(find(isnan(grads), 1)))
   keyboard; 
end