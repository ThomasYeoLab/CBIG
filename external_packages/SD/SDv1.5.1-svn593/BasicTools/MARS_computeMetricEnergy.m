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
function metric_energy = MARS_computeMetricEnergy(MARS_atlas, curr_sbjWarp, numWarpVerts, maxNeighbors) 

% metric_energy = MARS_computeMetricEnergy(MARS_atlas, curr_sbjWarp, numWarpVerts, maxNeighbors, C) 
%
% Compute Metric Distortion Energy (ignoring the partition function Z)

% C = 1000;
% %d_R2 is maxNeighbors x N
% %disp('Compute Distance Square between vertex and neighbors');
% %tic
% d_R2 = MARS_computeVertexDistSq2Nbors(int32(maxNeighbors), int32(numWarpVerts), int32(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors), single(curr_sbjWarp.curr_vertices));
% %toc
% 
% d_o2 = MARS_atlas.vertexDistSq2Nbors;
% 
% dR2_minus_do2 = d_R2 - d_o2;
% A = C.^4;
% 
% dR2_minus_do2_sq = dR2_minus_do2.^2;
% dR2_minus_do2_sq = (dR2_minus_do2_sq <= C^4).*dR2_minus_do2_sq + (dR2_minus_do2_sq > C^4).* (A*log(dR2_minus_do2_sq+single(eps)) - A*log(C^4) + C^4); 
% 
% % Note that for the entries corresponding to no neighbors, it's just 0, so
% % those won't affect anything
% %disp('Compute Metric Energy');
% %tic
% metric_energy = sum(sum(dR2_minus_do2_sq));
% %toc 


% This is a new metric energy which hopefully normalizes for multiscale.
num_neighbors_pairs = sum(sum(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors~=0));

d_R2 = MARS_computeVertexDistSq2Nbors(int32(maxNeighbors), int32(numWarpVerts), int32(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors), single(curr_sbjWarp.curr_vertices));
d_R = sqrt(d_R2);

d_o = sqrt(MARS_atlas.vertexDistSq2Nbors);

metric_energy = 983040/num_neighbors_pairs * sum(sum( ((d_R - d_o)./(d_o+eps)).^2)); %I added the 163842 so that it will be the same as the other energies...
