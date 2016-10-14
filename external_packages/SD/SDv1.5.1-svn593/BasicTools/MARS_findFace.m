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
function [face_indices, NViF] = MARS_findFace(points, mesh, seedVertices)
%function [face_indices, NViF] = MARS_findFace(points, mesh, seedVertices)
%returns: 
% face_index = The index of the face point falls into
% NV_ind = Index of nearest vertex in face
% distNV = Distance to nearest vertex in face
%MARS Project (Thomas Yeo and Mert Sabuncu), MIT, CSAIL, (c) 2006
%Author: Mert

if (nargin < 3)
    %disp('Finding seeds');
    [seedVertices, d] = MARS_findNVAux(single(points), single(mesh.vertices), int32(mesh.vertexNbors), int32(-1*ones(1, size(points,2))));
    %disp('Found seeds');
end
%disp('Find Face');
[face_indices, NViF] = MARS_findFaceAux(single(points), single(mesh.vertices), ...
    int32(mesh.faces), int32(mesh.vertexNbors), int32(mesh.vertexFaces), int32(seedVertices));
%disp('Found Face');
%return;
