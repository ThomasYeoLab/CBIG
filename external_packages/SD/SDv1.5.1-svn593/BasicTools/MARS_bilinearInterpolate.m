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
function vals = MARS_bilinearInterpolate(points, im)

% data_im = MARS_bilinearInterpolate(points, im)


radius = sqrt(sum(points.^2, 1));

offset_vec = zeros(1, size(points, 2));
offset_vec(abs(points(3, :) - 100) < 1e-5) = 1e-5;
offset_vec(abs(points(3, :) + 100) < 1e-5) = -1e-5;

theta_vec = mod(single(acos(points(3,:)./radius)) + offset_vec, single(pi));
phi_vec = mod(single(atan2(points(2,:), points(1,:))), single(2*pi));

num_data = int32(size(im, 1));
M = size(im, 2)-1;
N = size(im, 3)-1;

[phi, theta] = meshgrid(single(0):single(2*pi/N):single(2*pi), single(pi):single(-pi/M):single(0));
vals = zeros(num_data, length(theta_vec), 'single');
for i = 1:num_data
    temp = squeeze(im(i,:,:));
    vals(i,:) = MARS_interp2(phi, theta, temp, phi_vec, theta_vec, '*linear');
end


%vals = MARS_bilinearInterpolateAux(theta_vec, phi_vec, single(im), num_data, M, N, num_points);