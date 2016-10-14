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
function vec = SD_computeAtlas2SphereInvariantUpdate(current_data, warped_gauss_parms, parms)

% vec = SD_computeAtlas2SphereUpdate(current_data, warped_gauss_parms, parms)
%
% Written by Thomas, MIT

stdev = sqrt(warped_gauss_parms(2, :));
mesh = parms.meshes{parms.curr_level};


A_1 = (current_data - warped_gauss_parms(1, :))./sqrt(2)./stdev;

disp('Compute Weighted Gradient');
grad = MARS_linearInterpolateVertexAuxWGrad(mesh.vertices, mesh.faces, mesh.vertexFaces, mesh.faceAreas, warped_gauss_parms);
mu_grad = squeeze(grad(1, :, :));
stdev_grad = 0.5*squeeze(grad(2, :, :))./ repmat(stdev, 3, 1);

% Compute B_1
B_1 = - mu_grad./repmat((sqrt(2) * stdev), 3, 1) - repmat(A_1./stdev, 3, 1) .* stdev_grad;

% Compute B_2
B_2 = zeros(size(B_1));
greater_1_index = warped_gauss_parms(2, :) > 1+eps;
less_1_index = warped_gauss_parms(2, :) <= 1+eps;
B_2(:, greater_1_index) = 0.5*stdev_grad(:, greater_1_index)./repmat(sqrt(0.5*log(warped_gauss_parms(2, greater_1_index)).*warped_gauss_parms(2, greater_1_index)), 3, 1);

% Project onto tangent space.
b1 = [sum(B_1 .* mesh.e1, 1); sum(B_1 .* mesh.e2, 1)];
b2 = [sum(B_2 .* mesh.e1, 1); sum(B_2 .* mesh.e2, 1)];
b1b1T = OuterProduct2D(b1);
b2b2T = OuterProduct2D(b2);


% form the hessian
min_step = parms.meshes{parms.curr_level}.vertexDistSq2Nbors(1); %median(parms.meshes{parms.curr_level}.vertexDistSq2Nbors(parms.meshes{parms.curr_level}.vertexDistSq2Nbors~=0));
H = b1b1T + b2b2T + repmat(eye(2), [1 1 length(current_data)])./(parms.max_step^2 * min_step);

% Compute a2b2
a2b2 = stdev_grad/2./repmat(stdev, 3, 1);

% Comptue a1b1
a1b1 = repmat(A_1, 3, 1) .* B_1;

% Compute residual
res = a1b1 + a2b2; 

% Project residual
proj_res = [sum(res .* mesh.e1, 1); sum(res .* mesh.e2, 1)];


% Hessian Inverse
Hinv = -inverse2D(H);

% Calculate update 2 x N
vec2D = squeeze(sum(Hinv .* shiftdim(repmat(proj_res, [1 1 2]), 2), 2));
vec = mesh.e1 .* repmat(vec2D(1, :), 3, 1) + mesh.e2 .* repmat(vec2D(2, :), 3, 1);


function mat2D = OuterProduct2D(input_matrix)

% input  a 2 x N matrix
% output a 2 x 2 x N matrix, where 2x2 is obtained from an outer product
% from the input

mat2D = shiftdim(repmat(input_matrix, [1 1 2]), 2) .* shiftdim(repmat(input_matrix', [1 1 2]), 1);


function mat3D = OuterProduct(input_matrix)

% input  a 3 x N matrix
% output a 3 x 3 x N matrix, where 3x3 is obtained from an outer product
% from the input

mat3D = shiftdim(repmat(input_matrix, [1 1 3]), 2) .* shiftdim(repmat(input_matrix', [1 1 3]), 1);













































