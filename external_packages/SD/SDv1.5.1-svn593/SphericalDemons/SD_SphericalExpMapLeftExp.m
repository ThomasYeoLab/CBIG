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
function [best_warped_vertices, lowest_error] = SD_SphericalExpMapLeftExp(mesh, sbjWarp, parms, sbjMesh, basic_atlas, minDist2Nbors)

% [best_warped_vertices, lowest_error] = SD_SphericalExpMapLeftExp(mesh, sbjWarp, parms, sbjMesh, minDist2Nbors)
%
% SphericalExpMapLeftExp exponentiates sbjWarp.vec into displacement field on the
% sphere, with mesh.vertices replaced by parms.sbjWarp.curr_vertices as
% initialization. 
%
% exp(sbjWarp.vec) \circ parms.sbjWarp.curr_vertices 
%

% Compute input radius
input_radius = parms.radius;

if(~parms.use_uniform_mesh_ss)
    mesh.vertices = parms.sbjWarp.curr_vertices;
end
    
if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level && ~parms.use_uniform_mesh_ss)
    disp('Creating framework for fast latlon interpolation for scaling and squaring');
    if(parms.latlon_tight)
        num_per_side = floor(sqrt(size(basic_atlas.vertices, 2)/163842)*512);
        [temp, points, NF, mesh.numRows, mesh.numCols] = MARS_convertMesh2Image(mesh, sbjMesh.current_data, num_per_side, num_per_side);
    else
        [temp, points, NF, mesh.numRows, mesh.numCols] = MARS_convertMesh2Image(mesh, sbjMesh.current_data, 512, 512);
    end
    [mesh.weights, mesh.vno] = MARS_computeWeightsForInterpolation(mesh, points, NF);
end

% First compute smallest step.
if(nargin < 6)
   smallest_step = 1/(2^3); 
   max_velocity = sqrt(max(sum(sbjWarp.vec.^2, 1)));
   num_composition = ceil(log2(max_velocity/smallest_step));
else
   vec_stats = sqrt(sum(sbjWarp.vec.^2, 1)); 
   
   %we want vec / 2^n < d0/2^3
   num_composition = max(ceil(log2(max(vec_stats./minDist2Nbors)) + 3), 0);
   disp(['num_composition: ' num2str(num_composition)]);
end

%initialize output
lowest_error = parms.lowest_error;
best_warped_vertices = parms.sbjWarp.curr_vertices;

% set the first "first order integration"
warped_vertices = MARS_warpPointbyGradient(mesh.vertices, sbjWarp.vec/(2^num_composition));
    
% check if improvement is made.
if(parms.line_search_bool)
    
    if(~parms.use_uniform_mesh_ss)
        [test_warped_vertices, test_error] = computeObjFnDuringLeftExp(warped_vertices, parms, sbjMesh.current_data, basic_atlas);
    else
        temp_parms = parms;
        temp_parms.sbjWarp.curr_vertices = warped_vertices;
        temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
        temp_warped_vertices = temp_parms.sbjWarp.curr_vertices;
        
        temp_warped_vertices = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, mesh, temp_warped_vertices);
        [test_warped_vertices, test_error] = computeObjFnDuringLeftExp(temp_warped_vertices, parms, sbjMesh.current_data, basic_atlas);
    end

    lowest_error_index = -1;
    if(test_error < lowest_error)
        lowest_error_index = 0;
        lowest_error = test_error;
        best_warped_vertices = test_warped_vertices;
    end
end


uniform_mesh7 = mesh;
uniform_mesh7.vertices = warped_vertices;

% compute exp
for i = 1:(num_composition)
    if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level)% && ~parms.left_exp)
        warped_vertices = MARS_bilinearInterpolateDynamicImage(warped_vertices, mesh, warped_vertices);
    else
        warped_vertices = MARS_linearInterpolate(warped_vertices, mesh, warped_vertices);
    end
        
    new_radius = sqrt(sum(warped_vertices.^2, 1));
    warped_vertices = warped_vertices./repmat(new_radius, 3, 1) * input_radius;
    
    % check if improvement is made.
    if(parms.line_search_bool)
        if(~parms.use_uniform_mesh_ss)
            [test_warped_vertices, test_error] = computeObjFnDuringLeftExp(warped_vertices, parms, sbjMesh.current_data, basic_atlas);
        else
            temp_parms = parms;
            temp_parms.sbjWarp.curr_vertices = warped_vertices;
            temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
            temp_warped_vertices = temp_parms.sbjWarp.curr_vertices;
            
            temp_warped_vertices = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, mesh, temp_warped_vertices);
            [test_warped_vertices, test_error] = computeObjFnDuringLeftExp(temp_warped_vertices, parms, sbjMesh.current_data, basic_atlas);
        end

        if(test_error < lowest_error)
            lowest_error_index = i;
            lowest_error = test_error;
            best_warped_vertices = test_warped_vertices;
        end
    end
end

if(parms.line_search_bool)
    disp(num2str(lowest_error_index));
else
    if(~parms.use_uniform_mesh_ss)
        [best_warped_vertices, lowest_error] = computeObjFnDuringLeftExp(warped_vertices, parms, sbjMesh.current_data, basic_atlas);
    else

        temp_parms = parms;
        temp_parms.sbjWarp.curr_vertices = warped_vertices;
        temp_parms.sbjWarp = SD_smoothDeformationField(temp_parms, 1);
        temp_warped_vertices = temp_parms.sbjWarp.curr_vertices;

        temp_warped_vertices = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, mesh, temp_warped_vertices);
        [best_warped_vertices, lowest_error] = computeObjFnDuringLeftExp(temp_warped_vertices, parms, sbjMesh.current_data, basic_atlas);
    end
end


function [warped_vertices, error] = computeObjFnDuringLeftExp(warped_vertices, parms, current_data, basic_atlas)

if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level && ~parms.left_exp)
    warped_gauss_parms = MARS_bilinearInterpolate(warped_vertices, basic_atlas.im);
else
    warped_gauss_parms = MARS_linearInterpolate(warped_vertices, basic_atlas, [basic_atlas.mean; basic_atlas.variance]);
end
    
error = ComputeObjFnGivenData(current_data, warped_gauss_parms(1,:), warped_gauss_parms(2,:));






