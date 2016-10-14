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
function upsampled_warps = MARS_upsampleWarps(upsampledMesh, curr_warps)

% upsampledMesh is the mesh whose warps are to be upsampled
% curr_warps correspond to the warps of the downsampled version of
% upsampledMesh.
% Assume curr_warps = 3 x numPrevWarps
% Assume upsampled_warps(:, 1:numPrevWarps) = curr_warps;

numPrevWarps = size(curr_warps, 2);
numWarpVerts = size(upsampledMesh.vertices, 2);

upsampled_warps = zeros(3, numWarpVerts, 'single');
upsampled_warps(:, 1:numPrevWarps) = curr_warps;

vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(upsampledMesh.vertexNbors, 1)), int32(size(upsampledMesh.vertices, 2)), int32(upsampledMesh.vertexNbors), single(upsampledMesh.vertices));

%Only get those vertices we are interested in
vertexDistSq2Nbors = vertexDistSq2Nbors(:, numPrevWarps+1:end);
vertexNbors = upsampledMesh.vertexNbors(:, numPrevWarps+1:end);

% Set all the neighbors not from the previous resolution to 0.
vertexNbors(vertexNbors > numPrevWarps) = 0; 
vertexDistSq2Nbors = vertexDistSq2Nbors.*single(vertexNbors~=0);

% reshape vertexNbors to be of size 2 x (numWarpVerts - numPrevWarps)
tempVertexNbors = reshape(vertexNbors, numel(vertexNbors), 1);
tempVertexNbors = tempVertexNbors(tempVertexNbors~=0);
vertexNbors = reshape(tempVertexNbors, 2, numWarpVerts - numPrevWarps); %note that only 2 of the neighbors are from previous mesh, so this will work

% reshape vertexDistSq2Nbors to be of size 2 x (numWarpVerts - numPrevWarps) 
tempVertexDistSq2Nbors = reshape(vertexDistSq2Nbors, numel(vertexDistSq2Nbors), 1);
tempVertexDistSq2Nbors = tempVertexDistSq2Nbors(tempVertexDistSq2Nbors~=0);
vertexDistSq2Nbors = reshape(tempVertexDistSq2Nbors, 2, numWarpVerts - numPrevWarps); %note that only 2 of the neighbors are from previous mesh, so this will work
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);

interpolation_matrix = vertexDist2Nbors ./ repmat(sum(vertexDist2Nbors, 1), 2, 1);

warp_dim1_1 = curr_warps(1, vertexNbors(1,:));
warp_dim1_2 = curr_warps(1, vertexNbors(2,:));
avg_warp_dim1 = sum([warp_dim1_1; warp_dim1_2] .* interpolation_matrix,1);
upsampled_warps(1, numPrevWarps+1:end) = avg_warp_dim1;

warp_dim2_1 = curr_warps(2, vertexNbors(1,:));
warp_dim2_2 = curr_warps(2, vertexNbors(2,:));
avg_warp_dim2 = sum([warp_dim2_1; warp_dim2_2] .* interpolation_matrix,1);
upsampled_warps(2, numPrevWarps+1:end) = avg_warp_dim2;

warp_dim3_1 = curr_warps(3, vertexNbors(1,:));
warp_dim3_2 = curr_warps(3, vertexNbors(2,:));
avg_warp_dim3 = sum([warp_dim3_1; warp_dim3_2] .* interpolation_matrix,1);
upsampled_warps(3, numPrevWarps+1:end) = avg_warp_dim3;

input_radius = 100;
new_radius = sqrt(sum(upsampled_warps.^2, 1));
upsampled_warps = upsampled_warps./repmat(new_radius, 3, 1) * input_radius;
