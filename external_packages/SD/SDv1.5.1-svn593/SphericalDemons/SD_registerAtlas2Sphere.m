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
function [sbjWarp, harmonic_energy, lowest_error] = SD_registerAtlas2Sphere(sbjMesh, basic_atlas, parms)

% sbjWarp = SD_registerAtlas2Sphere(sbjMesh, basic_atlas, parms)
% 
% written by Thomas Yeo, MIT

input_radius = parms.radius;

mesh = parms.meshes{parms.curr_level};
origfaceAreas = MARS_computeMeshFaceAreas(int32(size(mesh.faces, 2)), int32(mesh.faces), single(mesh.vertices));


vertexDistSq2Nbors = mesh.vertexDistSq2Nbors;
vertexDistSq2Nbors(vertexDistSq2Nbors == 0) = max(vertexDistSq2Nbors(:));
minDist2Nbors = sqrt(min(vertexDistSq2Nbors, [], 1));

% Create Lat-lon interpolation
if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level) %&& (~parms.left_exp || (parms.left_exp && parms.use_uniform_mesh_ss)))
    disp('Creating framework for fast latlon interpolation');
    if(parms.latlon_tight)
        num_per_side = floor(sqrt(size(basic_atlas.vertices, 2)/163842)*512);
        basic_atlas.im = MARS_convertMesh2Image(basic_atlas, [basic_atlas.mean; basic_atlas.variance], num_per_side, num_per_side);
        [temp, points, NF, mesh.numRows, mesh.numCols] = MARS_convertMesh2Image(mesh, sbjMesh.current_data, num_per_side, num_per_side);
    else
        basic_atlas.im = MARS_convertMesh2Image(basic_atlas, [basic_atlas.mean; basic_atlas.variance], 512, 512);
        [temp, points, NF, mesh.numRows, mesh.numCols] = MARS_convertMesh2Image(mesh, sbjMesh.current_data, 512, 512);
    end
    [mesh.weights, mesh.vno] = MARS_computeWeightsForInterpolation(mesh, points, NF);
end
    
if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level)% && ~parms.left_exp)
    warped_gauss_parms = MARS_bilinearInterpolate(parms.sbjWarp.curr_vertices, basic_atlas.im);
else
    warped_gauss_parms = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, basic_atlas, [basic_atlas.mean; basic_atlas.variance]);
end

if(parms.multidim)
    BIG_NO = -1e6;
    index = find(mean(sbjMesh.current_data, 2) < 0.1*BIG_NO); %these structures do not exists
    sbjMesh.current_data(index, :) = warped_gauss_parms(index, :);
end
parms.lowest_error = ComputeObjFnGivenData(sbjMesh.current_data, warped_gauss_parms(1:end/2,:), warped_gauss_parms(end/2+1:end,:));
disp(['Before current scale begins, energy: ' num2str(parms.lowest_error)]);

for iter = 1:parms.iter(parms.curr_level)
    
    % Compute Update
    sbjWarp.vec = SD_computeAtlas2SphereInvariantUpdateRKHS2(sbjMesh.current_data, warped_gauss_parms, parms);
    
    sbjWarp.vec = MARS_projectGradOntoTangentPlane(mesh.vertices, sbjWarp.vec);

    if(sum(isnan(sbjWarp.vec(:))) > 0)
        keyboard;
    end
    
    % Smooth velocity updates
    if(parms.smooth_velocity)
        temp_parms = parms;
        temp_parms.sbjWarp = sbjWarp;
        sbjWarp = SD_smoothDeformationField(temp_parms, 1);
    end
    
    vec_stats = sqrt(sum(sbjWarp.vec.^2, 1));
    [parms.sbjWarp.curr_vertices, error_bs] = SD_SphericalExpMap(mesh, sbjWarp, parms, sbjMesh, basic_atlas, minDist2Nbors);
        
    new_radius = sqrt(sum(parms.sbjWarp.curr_vertices.^2, 1));
    parms.sbjWarp.curr_vertices = parms.sbjWarp.curr_vertices./repmat(new_radius, 3, 1) * input_radius;
     
    warp_bef = sqrt(sum((parms.sbjWarp.curr_vertices - mesh.vertices).^2, 1));
    
    % Compute if there's folding
    [folding_energy] = MARS_computeFoldingEnergyFast(parms.sbjWarp.curr_vertices, mesh);
    if(folding_energy > 0)
        if(parms.verbose)
            disp(['Folding energy (< 0 => folds exist): ' num2str(-1e6 * folding_energy)]);
        end
    end
            
    % Smooth with sigma
    parms.sbjWarp = SD_smoothDeformationField(parms);     
    
    % Unfold if necessary
    tempMesh = mesh;
    tempMesh.vertices = parms.sbjWarp.curr_vertices;
    tempMesh = MARS_simpleUnfoldMesh(tempMesh, 0);
    parms.sbjWarp.curr_vertices = tempMesh.vertices;
    
    if(~parms.minimal_verbosity || iter == parms.iter(parms.curr_level))
        % Compute area
        faceAreas = MARS_computeMeshFaceAreas(int32(size(mesh.faces, 2)), int32(mesh.faces), single(parms.sbjWarp.curr_vertices));
        maxFraction = max(faceAreas./origfaceAreas);
        minFraction = min(faceAreas./origfaceAreas);

        % Compute Energy after unfold
        warp = sqrt(sum((parms.sbjWarp.curr_vertices - mesh.vertices).^2, 1));

        disp(['max vec: ' num2str(max(vec_stats)) ', max warp bef smooth: ' num2str(max(warp_bef)) ', median warp bef smooth: ' num2str(median(warp_bef)) ...
            ', max warp after smooth: ' num2str(max(warp)) ', median warp after smooth: ' num2str(median(warp))]);

    end
    
    vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(mesh.vertexNbors, 1)), int32(size(parms.sbjWarp.curr_vertices, 2)), ...
                                    int32(mesh.vertexNbors), single(parms.sbjWarp.curr_vertices));                                
    harmonic_energy = sum((sqrt(vertexDistSq2Nbors(:)) - sqrt(mesh.vertexDistSq2Nbors(:))).^2)./sum(sum(mesh.vertexNbors~=0));
                                
    if(parms.latlon_interpolation && parms.curr_level >= parms.latlon_level) %&& ~parms.left_exp)
        warped_gauss_parms = MARS_bilinearInterpolate(parms.sbjWarp.curr_vertices, basic_atlas.im);
    else
        warped_gauss_parms = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, basic_atlas, [basic_atlas.mean; basic_atlas.variance]);
    end
    
    if(parms.multidim)
        sbjMesh.current_data(index, :) = warped_gauss_parms(index, :);
    end
    parms.lowest_error = ComputeObjFnGivenData(sbjMesh.current_data, warped_gauss_parms(1:end/2,:), warped_gauss_parms(end/2+1:end,:));
    
    if(~parms.minimal_verbosity || iter == parms.iter(parms.curr_level))
        disp([num2str(parms.curr_level) '.' num2str(iter) ': Harm Energy: ' num2str(harmonic_energy) ', Error (before smooth): ' num2str(error_bs) ...
            ', Error (after smooth): ' num2str(parms.lowest_error) ', maxArea: ' num2str(maxFraction) ', minArea: ' num2str(minFraction)]);

        disp('-----------------------------');
    end
    
end

sbjWarp = parms.sbjWarp;
lowest_error = parms.lowest_error;



