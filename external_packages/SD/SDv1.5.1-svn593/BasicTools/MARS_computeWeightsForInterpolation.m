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
function [weights, vno] = MARS_computeWeightsForInterpolation(mesh, points, NF)

% [weights, vno] = MARS_computeWeightsForInterpolation(mesh, points, NF)
%
% Written by Thomas Yeo
%
% This function is used by MARS_register2spheres. Aim is to compute weights
% so that we can quickly convert new data on the icosahedron onto lat-lon
% grid for fast interpolation.

weights = zeros(size(points), 'single');
vno = mesh.faces(:, NF);

v1 = mesh.vertices(:, mesh.faces(1, NF));
v2 = mesh.vertices(:, mesh.faces(2, NF));
v3 = mesh.vertices(:, mesh.faces(3, NF));

% First Compute normal to NF
v12 = v2 - v1;
v13 = v3 - v1;
normals = cross(v12, v13, 1);

% Now compute "projection". Note that this is the same as projectPoint in
% MARS_vec3D.h, which find intersection between vector representing a point
% and the triangle.
s = dot(v1, normals, 1)./dot(points, normals, 1);
proj_points = repmat(s, 3, 1) .* points;

% Compute 2*Areas
area1 = sqrt(sum(cross(proj_points - v2, proj_points - v3, 1).^2, 1));
area2 = sqrt(sum(cross(proj_points - v1, proj_points - v3, 1).^2, 1));
area3 = sqrt(sum(cross(proj_points - v1, proj_points - v2, 1).^2, 1));

total_area = area1+area2+area3;
weights(1, :) = area1./total_area;
weights(2, :) = area2./total_area;
weights(3, :) = area3./total_area;

