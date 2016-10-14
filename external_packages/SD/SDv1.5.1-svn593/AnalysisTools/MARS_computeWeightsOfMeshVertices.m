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
function weights = MARS_computeWeightsOfMeshVertices(MARS_sbjMesh, type)

% weights = MARS_computeWeightsOfMeshVertices(MARS_sbjMesh, type)
% if(strcmp(type, 'metric'))
%    vertices = MARS_sbjMesh.metricVerts;
% elseif(strcmp(type, 'spherical'))
%    vertices = MARS_sbjMesh.vertices;
% end

if(strcmp(type, 'metric'))
    vertices = MARS_sbjMesh.metricVerts;
elseif(strcmp(type, 'spherical'))
    vertices = MARS_sbjMesh.vertices;
end
    
vertexA = vertices(:, MARS_sbjMesh.faces(1,:));
vertexB = vertices(:, MARS_sbjMesh.faces(2,:));
vertexC = vertices(:, MARS_sbjMesh.faces(3,:));

vertexAB = vertexB - vertexA;
vertexAC = vertexC - vertexA;

areaFaces = cross(vertexAB, vertexAC, 1);
areaFaces = 0.5*sqrt(sum(areaFaces.^2, 1));
    
vertexFaces = MARS_sbjMesh.vertexFaces;
vertexFacesBool = double(vertexFaces ~= 0);
%face_count = sum(vertexFacesBool, 1);

vertexFaces(vertexFaces == 0) = 1; %avoiding indexing 0.

areaVertices = areaFaces(vertexFaces).*vertexFacesBool; % Multiplying will get rid of the empty faces

weights = 1/3 * sum(areaVertices, 1);

