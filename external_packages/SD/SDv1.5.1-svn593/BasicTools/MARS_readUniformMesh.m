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
function MARS_uniformMesh = MARS_readUniformMesh(mesh_dir, mesh_filename, radius)

% MARS_uniformMesh = MARS_readUniformMesh(mesh_dir, mesh_filename, radius)
%
%Given mesh dir and filename, read the uniform mesh and normalize it to
%radius

if(nargin < 3)
    radius = 100;
end

%read triangular mesh
[vertices, faces] = Read_Triangular_Mesh([mesh_dir '/' mesh_filename]);

%normalize radius
total_length = sqrt(sum(vertices.^2, 2));
vertices = vertices./repmat(total_length, 1, 3) * radius;

faces = int32(faces);

%Find vertexFaces
vertexFaces =  MARS_convertFaces2FacesOfVert(faces, int32(size(vertices, 1)));
num_per_vertex = length(vertexFaces)/size(vertices,1);
vertexFaces = reshape(vertexFaces, size(vertices,1), num_per_vertex);

%Find vertexNbors
vertexNbors = MARS_convertFaces2VertNbors(faces, int32(size(vertices,1)));
num_per_vertex = length(vertexNbors)/size(vertices,1);
vertexNbors = reshape(vertexNbors, size(vertices,1), num_per_vertex);

%Find vertexDistSq2Nbors
vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(vertexNbors', 1)), int32(size(vertices', 2)), int32(vertexNbors'), single(vertices'));

% Compute Face Areas.
faceAreas = MARS_computeMeshFaceAreas(int32(size(faces, 1)), int32(faces'), single(vertices'));  

%Create mesh structure and note the transpose!!
MARS_uniformMesh = struct('vertices', single(vertices)', 'faces', faces', 'vertexNbors', vertexNbors', 'vertexFaces', vertexFaces', 'vertexDistSq2Nbors', vertexDistSq2Nbors, 'faceAreas', faceAreas');

%Compute basis vectors
[MARS_uniformMesh.e1, MARS_uniformMesh.e2] = SD_findBasisVectors(MARS_uniformMesh.vertices);

% Originally planned to use below for fast gradient computation with matlab code. But
% current c implementation still faster than matlab code even though matlab
% code should be more efficient (less repeated computations...)

%Compute gradient vectors for each vertex
%MARS_uniformMesh.norm2base1 = ComputeNormal2TriBase(MARS_uniformMesh.vertices(:, faces(:, 1)), MARS_uniformMesh.vertices(:, faces(:, 2)), MARS_uniformMesh.vertices(:, faces(:, 3)));
%MARS_uniformMesh.norm2base2 = ComputeNormal2TriBase(MARS_uniformMesh.vertices(:, faces(:, 2)), MARS_uniformMesh.vertices(:, faces(:, 1)), MARS_uniformMesh.vertices(:, faces(:, 3)));
%MARS_uniformMesh.norm2base3 = ComputeNormal2TriBase(MARS_uniformMesh.vertices(:, faces(:, 3)), MARS_uniformMesh.vertices(:, faces(:, 2)), MARS_uniformMesh.vertices(:, faces(:, 1)));

%Compute projection matrix I - vn^T/<v,n> for each face connecting to each vertex.
%projMatCell is of length 6. Each projMat is 3 x 3 x numVertices.
%projMatCell{6}(:, :, 1:12) are zeros matrices because for uniform meshes,
%the first 12 vertices have only 5 faces. 
%MARS_uniformMesh.projMatCell = ComputeMatFromChainRuleFromProjOntoTriVertex(MARS_uniformMesh);









