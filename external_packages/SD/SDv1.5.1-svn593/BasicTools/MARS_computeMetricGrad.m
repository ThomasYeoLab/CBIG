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
function metric_grad = MARS_computeMetricGrad(MARS_atlas, curr_sbjWarp, numWarpVerts, maxNeighbors)

% metric_grad = MARS_computeMetricGrad(MARS_atlas, curr_sbjWarp, numWarpVerts, maxNeighbors) 
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
% % diffDataVertex2Nbors is 3 x maxNeighbors x numWarpVerts
% %disp('Compute Difference in Coordinates between vertex and neighbors: vertex - neighbor');
% %tic
% diffDataVertex2Nbors = MARS_computeDiffDataVertex2Nbors(int32(maxNeighbors), int32(numWarpVerts), int32(3), int32(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors), single(curr_sbjWarp.curr_vertices));
% %toc
% 
% dR2_minus_do2 = d_R2 - d_o2; 
% A = C.^4;
% 
% for i = 1:3
%    diffDataVertex2Nbors(i, :, :) = single(squeeze(diffDataVertex2Nbors(i, :, :)) .* ( dR2_minus_do2 .* (abs(dR2_minus_do2) <= C^2) ...
%                                                                                         + A.*(abs(dR2_minus_do2) > C^2) ./ (dR2_minus_do2 + single(eps))));    
% end
% 
% 
% 
% % Note that for the entries corresponding to no neighbors, it's just 0, so
% % those won't affect anything
% %disp('Compute the sum over neighbors: (d_R2 - d_o2) * (p_k - p_j)');
% %tic
% 
% metric_grad = squeeze(sum(diffDataVertex2Nbors, 2)); %now of size 3 x numWarpVerts
% %toc


% This is a new metric grad which hopefully normalizes for multiscale.
num_neighbors_pairs = sum(sum(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors~=0));

d_R2 = MARS_computeVertexDistSq2Nbors(int32(maxNeighbors), int32(numWarpVerts), int32(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors), single(curr_sbjWarp.curr_vertices));
d_R = sqrt(d_R2);

d_o = sqrt(MARS_atlas.vertexDistSq2Nbors);

% diffDataVertex2Nbors is 3 x maxNeighbors x numWarpVerts
%disp('Compute Difference in Coordinates between vertex and neighbors: vertex - neighbor');
diffDataVertex2Nbors = MARS_computeDiffDataVertex2Nbors(int32(maxNeighbors), int32(numWarpVerts), int32(3), int32(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors), single(curr_sbjWarp.curr_vertices));

% constant_factor is maxNeighbors x numWarpVerts
constant_factor = 2*(d_R - d_o)./(d_o+eps)./(d_R+eps);
for i = 1:3
   diffDataVertex2Nbors(i, :, :) = single(squeeze(diffDataVertex2Nbors(i, :, :))) .* constant_factor;    
end

metric_grad2 = 983040/num_neighbors_pairs * squeeze(sum(diffDataVertex2Nbors, 2)); %now of size 3 x numWarpVerts

% Now adding the adjustment factor due to the way we warp our vertices and
% then project onto the sphere

warp_factor = zeros(9, numWarpVerts, 'single'); % 9 x N
warp_factor(1,:) = curr_sbjWarp.curr_vertices(1,:) .* curr_sbjWarp.curr_vertices(1,:);
warp_factor(2,:) = curr_sbjWarp.curr_vertices(1,:) .* curr_sbjWarp.curr_vertices(2,:);
warp_factor(3,:) = curr_sbjWarp.curr_vertices(1,:) .* curr_sbjWarp.curr_vertices(3,:);
warp_factor(4,:) = warp_factor(2,:);
warp_factor(5,:) = curr_sbjWarp.curr_vertices(2,:) .* curr_sbjWarp.curr_vertices(2,:);
warp_factor(6,:) = curr_sbjWarp.curr_vertices(2,:) .* curr_sbjWarp.curr_vertices(3,:);
warp_factor(7,:) = warp_factor(3,:);
warp_factor(8,:) = warp_factor(6,:);
warp_factor(9,:) = curr_sbjWarp.curr_vertices(3,:) .* curr_sbjWarp.curr_vertices(3,:);

% I - 1/100^2 * x * x^T
warp_factor = repmat([1;0;0;0;1;0;0;0;1], 1, numWarpVerts) - warp_factor * 1/(100^2);

metric_grad = zeros(3, numWarpVerts, 'single');
metric_grad(1, :) = sum(warp_factor(1:3, :) .* metric_grad2, 1); 
metric_grad(2, :) = sum(warp_factor(4:6, :) .* metric_grad2, 1);
metric_grad(3, :) = sum(warp_factor(7:9, :) .* metric_grad2, 1);



% METRIC_CONST = single(20);
% RELATIVE_WEIGHT = 1e3;
% 
% temp = exp( min(METRIC_CONST * d_R2./(d_o2+eps), 50));
% 
% inf_index = find(temp == inf);
% if(~isempty(inf_index))
%    error('There is infinity!!');
% end
% 
% temp_index = (d_o2 == 0);
% mult_factor = (temp./(1+temp) - 1)./(d_o2+eps);
% mult_factor(temp_index) = 0;
% 
% diffDataVertex2Nbors = MARS_computeDiffDataVertex2Nbors(int32(maxNeighbors), int32(numWarpVerts), int32(3), int32(MARS_atlas.MARS_uniformMesh_SpatialPrior.vertexNbors), single(curr_sbjWarp.curr_vertices));
% 
% for i = 1:3
%     diffDataVertex2Nbors(i,:,:) = single(squeeze(diffDataVertex2Nbors(i, :, :))) .* mult_factor;
% end
% 
% augmented_metric_grad = RELATIVE_WEIGHT*squeeze(sum(diffDataVertex2Nbors, 2));
% 
% 
% metric_grad = metric_grad + augmented_metric_grad;