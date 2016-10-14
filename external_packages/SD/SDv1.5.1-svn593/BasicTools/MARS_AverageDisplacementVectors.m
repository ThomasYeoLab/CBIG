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
function smoothed_vertices = MARS_AverageDisplacementVectors(mesh, curr_vertices, var, num_iter)

% smoothed_vertices = MARS_AverageDisplacementVectors(uniform_mesh, curr_vertices, var, num_iter)
%
% Written by Thomas Yeo, MIT

%smoothed_vertices = MARS_findTangentVecPt1toPt2(mesh.vertices, curr_vertices);
smoothed_vertices = SD_TangentVecPt1toPt2Sine(mesh.vertices, curr_vertices);

for i = 1:num_iter
    smoothed_vertices = MARS_simpleAverageTangentVectors(mesh, smoothed_vertices, var);
    smoothed_vertices = MARS_projectGradOntoTangentPlane(mesh.vertices, smoothed_vertices);
end

%smoothed_vertices = MARS_warpPointbyTangentVec(mesh.vertices, smoothed_vertices);
smoothed_vertices = SD_warpPointbyTangentVecSine(mesh.vertices, smoothed_vertices);

% input_radius = 100;
% curr_vertices_sc =  [acos(curr_vertices(3, :)./input_radius);  atan2(curr_vertices(2,:), curr_vertices(1,:))]; % [theta; phi]
% curr_vertices_sc(1, :) = mod(curr_vertices_sc(1, :) + 2*pi, 2*pi);
% curr_vertices_sc(2, :) = mod(curr_vertices_sc(2, :) + 2*pi, 2*pi);
% 
% for i = 1:num_iter
%     curr_vertices_sc = MARS_simpleAverageSC(mesh, curr_vertices_sc);
% end
% 
% curr_vertices_sc(1, :) = mod(curr_vertices_sc(1, :) + 2*pi, 2*pi);
% curr_vertices_sc(2, :) = mod(warped_vertices_sc(2, :) + 2*pi, 2*pi);
% smoothed_vertices = [input_radius*cos(curr_vertices_sc(2,:)).*sin(curr_vertices_sc(1,:)); ...
%                         input_radius*sin(curr_vertices_sc(2,:)).*sin(curr_vertices_sc(1,:)); ...
%                         input_radius*cos(curr_vertices_sc(1,:))]; 