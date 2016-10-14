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
function vec = SD_computeAtlas2SphereInvariantUpdateRKHS2(current_data, warped_gauss_parms, parms)

% vec = SD_computeAtlas2SphereUpdate(current_data, warped_gauss_parms, parms)
%
% Written by Thomas, MIT

if(size(current_data, 1) == 1)
   parms.multidim = 0; 
end

stdev = sqrt(warped_gauss_parms(end/2+1:end, :));
mesh = parms.meshes{parms.curr_level};

A_1 = (current_data - warped_gauss_parms(1:end/2, :))./sqrt(2)./stdev;
A_1_3D = zeros(size(A_1, 1), 3, size(A_1, 2));
A_1_3D(:, 1, :) = A_1;
A_1_3D(:, 2, :) = A_1;
A_1_3D(:, 3, :) = A_1;
A_1_3D = squeeze(A_1_3D);

if(parms.verbose)
    disp('Compute Weighted Gradient');
end
grad = MARS_linearInterpolateVertexAuxWGrad(mesh.vertices, mesh.faces, mesh.vertexFaces, mesh.faceAreas, warped_gauss_parms);
mu_grad = squeeze(grad(1:end/2, :, :));

stdev3D = zeros(size(stdev, 1), 3, size(stdev, 2));
stdev3D(:, 1, :) = stdev;
stdev3D(:, 2, :) = stdev;
stdev3D(:, 3, :) = stdev;
stdev3D = squeeze(stdev3D);
stdev_grad = 0.5*squeeze(grad(end/2+1:end, :, :))./ stdev3D;

% Compute B_1
B_1 = - mu_grad./(sqrt(2) * stdev3D) - A_1_3D./stdev3D .* stdev_grad;

% Compute B_2
B_2 = stdev_grad./stdev3D;

% Compute dist grad
if(parms.rkhs)
    if(parms.geodesic)
        dist_grad = MARS_linearInterpolateVertexDisplacementWGrad(mesh.vertices, mesh.faces, mesh.vertexFaces, mesh.faceAreas, parms.sbjWarp.curr_vertices, int32(1));
    else
        dist_grad = MARS_linearInterpolateVertexDisplacementWGrad(mesh.vertices, mesh.faces, mesh.vertexFaces, mesh.faceAreas, parms.sbjWarp.curr_vertices, int32(0));
    end
end


% This is wrong... needs to smooth Hessian. Not just smooth derivatives...
% if(parms.smooth_velocity)
%     temp_parms = parms;
%     
%     temp_parms.sbjWarp.vec = B_1;
%     temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
%     B_1 = temp_parms.sbjWarp.vec;
%     
%     temp_parms.sbjWarp.vec = B_2;
%     temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
%     B_2 = temp_parms.sbjWarp.vec;
%     
%     if(parms.rkhs)
%         temp_parms.sbjWarp.vec = squeeze(dist_grad(1, :, :));
%         temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
%         dist_grad(1, :, :) = temp_parms.sbjWarp.vec;
% 
%         temp_parms.sbjWarp.vec = squeeze(dist_grad(2, :, :));
%         temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
%         dist_grad(2, :, :) = temp_parms.sbjWarp.vec;
% 
%         temp_parms.sbjWarp.vec = squeeze(dist_grad(3, :, :));
%         temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
%         dist_grad(3, :, :) = temp_parms.sbjWarp.vec;
%     end
% end

% Project B_1 onto tangent space and estimate hessian due to B1
e1 = shiftdim(repmat(mesh.e1, [1 1 size(A_1, 1)]), 2);
e2 = shiftdim(repmat(mesh.e2, [1 1 size(A_1, 1)]), 2);

b1 = zeros(size(current_data,1), 2, size(current_data, 2));
if(parms.multidim)
    b1(:, 1, :) = sum(B_1 .* e1, 2);
    b1(:, 2, :) = sum(B_1 .* e2, 2);
else
    b1(:, 1, :) = sum(B_1 .* e1, 1);
    b1(:, 2, :) = sum(B_1 .* e2, 1);
end
b1b1T = OuterProduct2D_2(b1);

% Project B_2 onto tangent space and estimate Hessian due to B1
b2 = zeros(size(current_data,1), 2, size(current_data, 2));
if(parms.multidim)
    b2(:, 1, :) = sum(B_2 .* e1, 2);
    b2(:, 2, :) = sum(B_2 .* e2, 2);
else
    b2(:, 1, :) = sum(B_2 .* e1, 1);
    b2(:, 2, :) = sum(B_2 .* e2, 1);
end

if(parms.approxH)
    approxH = -OuterProduct2D_2(b2);
else
    approxH = 0;
end

% Estimate Hessian due to dist_grad
if(parms.rkhs)
    proj_dist_grad = zeros(3, 2, length(stdev));

    temp_grad1 = dist_grad;
    temp_grad1(1, :, :) = squeeze(temp_grad1(1, :, :)) .* mesh.e1;
    temp_grad1(2, :, :) = squeeze(temp_grad1(2, :, :)) .* mesh.e1;
    temp_grad1(3, :, :) = squeeze(temp_grad1(3, :, :)) .* mesh.e1;

    temp_grad2 = dist_grad;
    temp_grad2(1, :, :) = squeeze(temp_grad2(1, :, :)) .* mesh.e2;
    temp_grad2(2, :, :) = squeeze(temp_grad2(2, :, :)) .* mesh.e2;
    temp_grad2(3, :, :) = squeeze(temp_grad2(3, :, :)) .* mesh.e2;

    proj_dist_grad(:, 1, :) = squeeze(sum(temp_grad1, 2));
    proj_dist_grad(:, 2, :) = squeeze(sum(temp_grad2, 2));

    proj_hess = zeros(2, 2, length(stdev));
    proj_hess(1, 1, :) = sum(squeeze(proj_dist_grad(:, 1, :)) .* squeeze(proj_dist_grad(:, 1, :)), 1);
    proj_hess(2, 2, :) = sum(squeeze(proj_dist_grad(:, 2, :)) .* squeeze(proj_dist_grad(:, 2, :)), 1);
    proj_hess(1, 2, :) = sum(squeeze(proj_dist_grad(:, 1, :)) .* squeeze(proj_dist_grad(:, 2, :)), 1);
    proj_hess(2, 1, :) = proj_hess(1,2,:);
    proj_hess = proj_hess / 512 * ((length(stdev)/2562)^2);
end

% form the hessian
min_step = parms.meshes{parms.curr_level}.vertexDistSq2Nbors(1); %median(parms.meshes{parms.curr_level}.vertexDistSq2Nbors(parms.meshes{parms.curr_level}.vertexDistSq2Nbors~=0));

id = repmat(eye(2), [1 1 length(current_data)]);
if(parms.rkhs)
    H = b1b1T + 0.5*(-approxH) + 0.5*(proj_hess + id)./(parms.max_step^2 * min_step);
else
    H = b1b1T + 0.5*(-approxH) + id./(parms.max_step^2 * min_step);
end

% Comptue a1b1
% Compute residual
if(parms.multidim)
    a1b1 = squeeze(sum(A_1_3D .* B_1, 1));
    res = a1b1 + 0.5*squeeze(sum(B_2, 1)); 
else
    a1b1 = A_1_3D .* B_1;
    res = a1b1 + 0.5*B_2; 
end
    
% Project residual
proj_res = [sum(res .* mesh.e1, 1); sum(res .* mesh.e2, 1)];


% Hessian Inverse
Hinv = -inverse2D(H);

% Calculate update 2 x N
vec2D = squeeze(sum(Hinv .* shiftdim(repmat(proj_res, [1 1 2]), 2), 2));
vec = mesh.e1 .* repmat(vec2D(1, :), 3, 1) + mesh.e2 .* repmat(vec2D(2, :), 3, 1);


function mat2D = OuterProduct2D_2(input_matrix)

% input  a D x 2 x N matrix
% output a 2 x 2 x N matrix, where 2x2 is obtained from an outer product
% from the input

mat2D = zeros(2, 2, size(input_matrix, 3));
mat2D(1, 1, :) = squeeze(sum(input_matrix(:, 1, :) .* input_matrix(:, 1, :), 1));
mat2D(2, 2, :) = squeeze(sum(input_matrix(:, 2, :) .* input_matrix(:, 2, :), 1));
mat2D(1, 2, :) = squeeze(sum(input_matrix(:, 1, :) .* input_matrix(:, 2, :), 1));
mat2D(2, 1, :) = mat2D(1, 2, :);


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













































