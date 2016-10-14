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
function image_grad = SD_computeGradAtVerticesWithBias(mesh, data, pos_bool)

% image_grad = SD_computeGradAtVerticesWithBias(mesh, data, pos_bool)
%
% Written by Thomas Yeo
% Assume data = 1 x N right now. 
%
% Assume pos_bool consists of ones and minus ones. ones implies that we are
% looking for biggest values.

vertexNbors = mesh.vertexNbors;
zero_index = (vertexNbors == 0);
vertexNbors(zero_index) = 1;

neighbors_data = data(vertexNbors);

% switch signs of neighbors which you want to find smallest values.
vertexDistSq2Nbors = mesh.vertexDistSq2Nbors;
vertexDistSq2Nbors(vertexDistSq2Nbors == 0) = max(vertexDistSq2Nbors(:));
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);


neighbors_gradient = (neighbors_data - repmat(data, size(vertexNbors, 1), 1))./vertexDist2Nbors; 
neighbors_gradient_temp = neighbors_gradient .* repmat(pos_bool, size(vertexNbors, 1), 1);
neighbors_gradient_temp(zero_index) = -inf;


% First find "largest" neighbors
[max_data, max_index] = max(neighbors_gradient_temp, [], 1);
index_grad_neighbors = sub2ind(size(vertexNbors), max_index, 1:size(vertexNbors, 2));
max_vertex = vertexNbors(index_grad_neighbors);

grad_direction = mesh.vertices(:, max_vertex) - mesh.vertices;
grad_direction = grad_direction ./ repmat(sqrt(sum(grad_direction.^2, 1)), 3, 1); %unit gradient direction

image_grad = repmat(neighbors_gradient(index_grad_neighbors), 3, 1) .* grad_direction;

% So far, we have computed the image gradient as we move along the
% triangle. But the triangle does not lie on the sphere.
% So we need to take into account the projection.
%
% see my notebook for derivations
% A = I - 1/<v, n> vn^T 
% image_grad * A = A^T * image_grad' 
%
% Note that in this case, we are moving along the edge of two triangle. We
% can just use any normal that is perpendicular to p1 and it's neighbor.
%

% First first direction to neighbors;
% vec2Nbor = mesh.vertices(:, max_vertex) - mesh.vertices;
% NormVec = cross(mesh.vertices, vec2Nbor, 1); 
% NormVec = cross(NormVec, vec2Nbor, 1);
% 
% 
% nTv = reshape(repmat(NormVec, [3 1 1]), 3, 3, size(NormVec, 2)) .* shiftdim(repmat(mesh.vertices, [1 1 3]), 2);
% AT = repmat(eye(3, 3), [1 1 size(mesh.vertices, 2)]) - nTv./squeeze(shiftdim(repmat(dot(mesh.vertices, NormVec, 1), [3 1 3]), 2)); 
% 
% image_grad = squeeze(sum(AT .* shiftdim(repmat(image_grad, [1 1 3]), 2), 2));

image_grad = MARS_projectGradOntoTangentPlane(mesh.vertices, image_grad);










