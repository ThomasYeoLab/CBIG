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
function MARS_sbjMesh = MARS_simpleUnfoldMesh(MARS_sbjMesh, reverse_face_bool, max_iter)

% MARS_sbjMesh = MARS_simpleUnfoldMesh(MARS_sbjMesh, reverse_face_bool, max_iter)
% reverse_face_bool = 1 if orientation is opposite from that of icosahedron

if(nargin < 3)
   max_iter = 1e5; 
end

if(reverse_face_bool)
    faces = MARS_sbjMesh.faces;
    MARS_sbjMesh.faces = [faces(3,:); faces(2,:); faces(1,:)];
end


[folding_energy, folded_vertices_list] = MARS_computeFoldingEnergyFast(MARS_sbjMesh.vertices, MARS_sbjMesh);
folding_energy = -1e6 * folding_energy;
   
if(folding_energy < 0)  
    %tic
    disp('Unfolding Mesh!!');
    sbjWarp = struct('curr_vertices', MARS_sbjMesh.vertices);
    vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)), int32(size(MARS_sbjMesh.vertices, 2)), int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.vertices));
    MARS_atlas = struct('vertexDistSq2Nbors', single(vertexDistSq2Nbors), 'MARS_uniformMesh_SpatialPrior', MARS_sbjMesh);

    MARS_sbjMesh.vertices = MARS_unfoldMesh(MARS_atlas, sbjWarp, folding_energy, 0.3 * sqrt(min(MARS_atlas.vertexDistSq2Nbors(1,:))), 4, 1e6, folded_vertices_list, max_iter);
end

%switch it back!!    
if(reverse_face_bool)
    faces = MARS_sbjMesh.faces;
    MARS_sbjMesh.faces = [faces(3,:); faces(2,:); faces(1,:)];
end