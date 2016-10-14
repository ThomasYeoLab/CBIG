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
function [best_warped_vertices, lowest_error] = SD_SphericalExpMap(mesh, sbjWarp, parms, sbjMesh, basic_atlas, minDist2Nbors)

% [best_warped_vertices, lowest_error] = SD_SphericalExpMapLeftExp(mesh, sbjWarp, parms, sbjMesh, minDist2Nbors)
%
% SphericalExpMap exponentiates sbjWarp.vec into displacement field on the
% sphere, with mesh.vertices as initialization. It then composes this 
% displacement field with sbjWarp.curr_vertices.
%
% sbjWarp.curr_vertices \circ exp(sbjWarp.vec)
%
% Note that algorithm assumes inter-vertex distance is about 1.
% Therefore if step-size is less than 0.5, things should stay
% diffeomorphic. 


% Compute input radius
input_radius = parms.radius;

% First compute smallest step.
if(nargin < 6)
   smallest_step = 1/(2^3); 
   max_velocity = sqrt(max(sum(sbjWarp.vec.^2, 1)));
   num_composition = ceil(log2(max_velocity/smallest_step));
else
   vec_stats = sqrt(sum(sbjWarp.vec.^2, 1)); 
   
   %we want vec / 2^n < d0/2^3
   num_composition = max(ceil(log2(max(vec_stats./minDist2Nbors)) + 3), 0);
   if(parms.verbose)
       disp(['scaling and squaring: num_composition: ' num2str(num_composition)]);
   end
end

% set the first "first order integration" 
% Note that using MARS_warpPointbyGradient or MARS_warpPointbyTangentVec
% results in very minor differences. This is simply because the chart we use and the 
% exponential chart agrees up to the first order at the origin. Therefore,
% for small velocity (v divided 2^N), they are essentially the same. 
%
% For example, in the trial run we did, with the chart we use, we get: 
%       Orig Error: 88016.0781, Final Error: 14772.7529, Harm Energy: 0.051755
% 
% With the exponential chart, we get:
%       Orig Error: 88016.0781, Final Error: 14772.7666, Harm Energy: 0.051756
%
% Hence, they are the same up to 5 significant figures.
% The default is set to MARS_warpPointbyGradient since it is faster.
warped_vertices = MARS_warpPointbyGradient(mesh.vertices, sbjWarp.vec/(2^num_composition));
%warped_vertices = MARS_warpPointbyTangentVec(mesh.vertices, sbjWarp.vec/(2^num_composition));
    
uniform_mesh7 = mesh;
uniform_mesh7.vertices = warped_vertices;

% compute exp
for i = 1:(num_composition)
    
    if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level)
        warped_vertices = MARS_bilinearInterpolateDynamicImage(warped_vertices, mesh, warped_vertices);
    else
        warped_vertices = MARS_linearInterpolate(warped_vertices, mesh, warped_vertices);
    end
        
    new_radius = sqrt(sum(warped_vertices.^2, 1));
    warped_vertices = warped_vertices./repmat(new_radius, 3, 1) * input_radius;
    
end

[best_warped_vertices, lowest_error] = computeObjFnDuringExp(mesh, warped_vertices, parms, sbjMesh.current_data, basic_atlas);


function [warped_vertices, error] = computeObjFnDuringExp(mesh, warped_vertices, parms, current_data, basic_atlas)

if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level)
    warped_vertices = MARS_bilinearInterpolateDynamicImage(warped_vertices, mesh, parms.sbjWarp.curr_vertices);
else
    warped_vertices = MARS_linearInterpolate(warped_vertices, mesh, parms.sbjWarp.curr_vertices);
end

new_radius = sqrt(sum(warped_vertices.^2, 1));
warped_vertices = warped_vertices./repmat(new_radius, 3, 1) * parms.radius;


if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level)
    warped_gauss_parms = MARS_bilinearInterpolate(warped_vertices, basic_atlas.im);
else
    warped_gauss_parms = MARS_linearInterpolate(warped_vertices, basic_atlas, [basic_atlas.mean; basic_atlas.variance]);
end
    
error = ComputeObjFnGivenData(current_data, warped_gauss_parms(1:end/2,:), warped_gauss_parms(end/2+1:end,:));






