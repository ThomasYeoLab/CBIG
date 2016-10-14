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
function [dist_transform, LogOdds, normalized_prob_mat] = MARS_computeLogOdds(MARS_sbjMesh, alpha, BIG_NO)

% LogOdds = MARS_computeLogOdds(MARS_sbjMesh, alpha)
% Note that alpha scales the sign distance function and hence determines
% how much confidence we place in the boundaries.

num_structures = MARS_sbjMesh.MARS_ct.numEntries;
num_vertices = int32(size(MARS_sbjMesh.vertices, 2));
maxNeighbors = int32(size(MARS_sbjMesh.vertexNbors, 1));
dist_transform = zeros(num_structures, num_vertices);

vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(MARS_sbjMesh.vertexNbors, 1)), int32(size(MARS_sbjMesh.vertices, 2)), int32(MARS_sbjMesh.vertexNbors), single(MARS_sbjMesh.vertices));
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);

for i = 1:num_structures
    
    objVec = (MARS_sbjMesh.MARS_label == i); 
    
    if(~isempty(find(objVec == 1)))
        [boundaryVec] = ConvertObjVec2Boundary(MARS_sbjMesh, objVec);

        temp_dist_transform = MARS_DT_Boundary(int32(boundaryVec), num_vertices, maxNeighbors, MARS_sbjMesh.vertexNbors, double(vertexDist2Nbors)); %min_heap assumes double
        temp_dist_transform(~objVec) = -temp_dist_transform(~objVec);
        dist_transform(i, :) = temp_dist_transform;
    else
        dist_transform(i, :) = -BIG_NO;
    end
end

dist_transform = dist_transform*alpha;

%Convert To Probability
prob_mat = 1./(1 + exp(-dist_transform)) + eps;
total_prob = sum(prob_mat, 1);

normalized_prob_mat = prob_mat./repmat(total_prob, num_structures, 1);

LogOdds = log(normalized_prob_mat./(1-normalized_prob_mat));










