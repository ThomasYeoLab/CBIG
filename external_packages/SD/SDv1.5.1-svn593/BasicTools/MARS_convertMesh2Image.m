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
function [im, points, NF, numRows, numCols] = MARS_convertMesh2Image(mesh, data, M, N, radius)

% im = ConvertMesh2Image(mesh, data, M, N, radius)
% Convert data that lies on mesh into a (M+1:N+1) matrix 
% where the first and last column are the same.

if(nargin < 5)
   radius = 100; 
end

[phi, theta] = meshgrid(0:2*pi/N:2*pi, pi:-pi/M:0);

num_verts = (M+1)*(N+1);

phi_vec = reshape(phi, 1, num_verts);
theta_vec = reshape(theta, 1, num_verts);

points = single(radius*[cos(phi_vec).*sin(theta_vec); sin(phi_vec).*sin(theta_vec); cos(theta_vec)]);

[im, NV, NF] = MARS_linearInterpolate(points, mesh, data);

num_data = size(data, 1);

im = reshape(im, num_data, M+1, N+1);
numRows = M+1;
numCols = N+1;

