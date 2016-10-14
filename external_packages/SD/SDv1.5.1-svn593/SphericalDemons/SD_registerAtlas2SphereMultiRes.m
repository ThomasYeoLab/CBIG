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
function parms = SD_registerAtlas2SphereMultiRes(sbjMesh, SD_atlas, parms)

% parms = SD_registerAtlas2SphereMultiRes(sbjMesh, SD_atlas, parms)
%
% Diffeomorphic spherical registration of atlas to spherical subject mesh in a
% multi-level fashion. Atlas is warped.
% 
% written by Thomas Yeo, MIT
%
% Assumes that initialization from parms.sbjWarp and that parms.sbjWarp is
% the same SIZE as sbjMesh.

if(parms.verbose)
    disp(parms);
end
    
if(SD_atlas.multidim)
    BIG_NO = -1e6;
    index = find(mean(sbjMesh.data, 2) > 0.1*BIG_NO); %these structures do  exist
    sbjMesh.data = sbjMesh.data(index, :);
    for j = 1:length(SD_atlas.basic_atlas)
       SD_atlas.basic_atlas{j}.mean = SD_atlas.basic_atlas{j}.mean(index, :);
       SD_atlas.basic_atlas{j}.variance = SD_atlas.basic_atlas{j}.variance(index, :);
    end
end

% Initialize warp if necessary
if(isempty(parms.sbjWarp))
   no_initialization = 1; 
   sbjWarp.curr_vertices = sbjMesh.vertices;
   parms.sbjWarp  = sbjWarp;
else
   no_initialization = 0; 
end

if(size(parms.sbjWarp.curr_vertices, 2) ~= size(sbjMesh.vertices, 2))
        error('Dimension of initialized warp not equal to sbjMesh');
end

% Compute energy before registration
if(SD_atlas.multidim)
    sbjMesh.current_data = sbjMesh.data;
else
    sbjMesh.current_data = sbjMesh.data(end,:);
end
orig_objfn = ComputeObjFnGivenMeshes(sbjMesh, SD_atlas, parms);

RotationSearchWidth = parms.rotate_width;
numInterval = parms.rotate_interval;
for level = 1:length(parms.meshes)
    
    parms.curr_level = level;

    if(level == 1)
        % need to downsample warps
        disp('Downsampling sbjWarp');
        parms.sbjWarp.curr_vertices = MARS_linearInterpolate(parms.meshes{1}.vertices, sbjMesh, parms.sbjWarp.curr_vertices);
    else
        disp('Interpolating downsampled_sbjWarp to higher resolution');
        upsampledMesh = parms.meshes{level};
        parms.sbjWarp.curr_vertices = MARS_upsampleWarps(upsampledMesh, parms.sbjWarp.curr_vertices);
    end
        
    % Unfold if necessary
    tempMesh = parms.meshes{level};
    tempMesh.vertices = parms.sbjWarp.curr_vertices;
    tempMesh = MARS_simpleUnfoldMesh(tempMesh, 0);
    parms.sbjWarp.curr_vertices = tempMesh.vertices;
    
    % Get sbjMesh data of current resolution
    if(SD_atlas.multidim)
        sbjMesh.current_data = MARS_linearInterpolate(parms.meshes{level}.vertices, sbjMesh, sbjMesh.data);        
    else
        sbjMesh.current_data = MARS_linearInterpolate(parms.meshes{level}.vertices, sbjMesh, sbjMesh.data(level, :));        
    end
        
    if(parms.iter(parms.curr_level) > 0)   
        if(parms.movie_flag && ~SD_atlas.multidim)
            warped_mean = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, SD_atlas.basic_atlas{level}, SD_atlas.basic_atlas{level}.mean);
            temp_data = sbjMesh.current_data;

            max_value = max(max(temp_data), max(warped_mean));
            min_value = min(min(temp_data), min(warped_mean));
            
            temp_data2 = temp_data; temp_data2(1)=max_value; temp_data2(2) = min_value;
            figure; set(gcf, 'renderer', 'zbuffer'); TrisurfMeshData(parms.meshes{level}, temp_data2); shading interp; colorbar; view(250, 20); 
            title('Subject data'); saveas(gcf, ['temp' num2str(level*4 - 3, '%.3d') '.tiff']);
            
            warped_mean2 = warped_mean; warped_mean2(1)=max_value; warped_mean2(2) = min_value;
            figure; set(gcf, 'renderer', 'zbuffer'); TrisurfMeshData(parms.meshes{level}, warped_mean2); shading interp; colorbar; view(250, 20); 
            title(['Atlas mean at start of multiresolution level ' num2str(level)]); saveas(gcf, ['temp' num2str(level*4 - 2, '%.3d') '.tiff']);
        end

        parms.multidim = SD_atlas.multidim;
        
        % Perform Rotation Alignment
        if(parms.rotate)
            disp('Align by rotation');
            [tmp.sbjWarp, Center1, Center2, Center3, pe, ce] = SD_rotateAtlas2Sphere(sbjMesh, SD_atlas.basic_atlas{level}, parms, RotationSearchWidth, numInterval);
            disp(['Rotate by ' num2str(Center1) ', ' num2str(Center2) ', ' num2str(Center3) '( prev energy: ' num2str(pe) ', curr energy: ' num2str(ce) ')']);
            parms.sbjWarp.curr_vertices = tmp.sbjWarp.curr_vertices;
            %if(level == 1 && no_initialization) % if no warp at first, better to just adjust original warps.
            %    sbjMesh.vertices = MARS_zrotate(MARS_yrotate(MARS_zrotate(sbjMesh.vertices, Center3), Center2), Center1);
            %    if(SD_atlas.multidim)
            %        sbjMesh.current_data = MARS_linearInterpolate(parms.meshes{level}.vertices, sbjMesh, sbjMesh.data);
            %    else
            %        sbjMesh.current_data = MARS_linearInterpolate(parms.meshes{level}.vertices, sbjMesh, sbjMesh.data(level, :));
            %    end
            %else %if there's already warps, should just update warp vertices.
            %    parms.sbjWarp.curr_vertices = tmp.sbjWarp.curr_vertices;
            %end

            if(parms.rotate_multiscale)
                RotationSearchWidth = RotationSearchWidth/2;
                numInterval = max(numInterval/2, 2);
            end
            
            if(parms.movie_flag  && ~SD_atlas.multidim)
                warped_mean = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, SD_atlas.basic_atlas{level}, SD_atlas.basic_atlas{level}.mean);
                
                warped_mean2 = warped_mean; warped_mean2(1)=max_value; warped_mean2(2) = min_value;
                figure; set(gcf, 'renderer', 'zbuffer'); TrisurfMeshData(parms.meshes{level}, warped_mean2); shading interp; colorbar; view(250, 20); 
                title(['Atlas mean after rotation ' num2str(level)]); saveas(gcf, ['temp' num2str(level*4 - 1, '%.3d') '.tiff']);
            end
        end
            
        % Perform Single Res Registration
        disp(['================== Perform Single Res Registration: ' num2str(level) ' ==================']);
        [parms.sbjWarp, parms.harmonic_energy_bef, parms.error_bef_interp_warps] = SD_registerAtlas2Sphere(sbjMesh, SD_atlas.basic_atlas{level}, parms);

        if(parms.movie_flag  && ~SD_atlas.multidim)
            warped_mean = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, SD_atlas.basic_atlas{level}, SD_atlas.basic_atlas{level}.mean);
            
            warped_mean2 = warped_mean; warped_mean2(1)=max_value; warped_mean2(2) = min_value;
            figure; set(gcf, 'renderer', 'zbuffer'); TrisurfMeshData(parms.meshes{level}, warped_mean2); shading interp; colorbar; view(250, 20); 
            title(['Atlas mean at end of multiresolution level ' num2str(level)]); saveas(gcf, ['temp' num2str(level*4, '%.3d') '.tiff']);
            close all;
        end  
    end
end


% temp_mesh = parms.meshes{parms.curr_level};
% temp_mesh.vertices = parms.sbjWarp.curr_vertices;
% parms.sbjWarp.curr_vertices = MARS_linearInterpolate(SD_atlas.basic_atlas{level}.vertices, temp_mesh, parms.meshes{level}.vertices);
% 
% vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(temp_mesh.vertexNbors, 1)), int32(size(parms.sbjWarp.curr_vertices, 2)), ...
%      int32(temp_mesh.vertexNbors), single(parms.sbjWarp.curr_vertices));
% parms.harmonic_energy = sum((sqrt(vertexDistSq2Nbors(:)) - sqrt(temp_mesh.vertexDistSq2Nbors(:))).^2)./sum(sum(temp_mesh.vertexNbors~=0));
% 
% parms.lowest_error = ComputeSDObjFnGivenMeshes(sbjMesh, SD_atlas, parms);
% disp(['Orig Error: ' num2str(orig_objfn) ', Final Error: ' num2str(parms.lowest_error) ', Harm Energy: ' num2str(parms.harmonic_energy)]);
% 


parms.sbjWarp.curr_vertices = MARS_linearInterpolate(sbjMesh.vertices, parms.meshes{parms.curr_level}, parms.sbjWarp.curr_vertices);
new_radius = sqrt(sum(parms.sbjWarp.curr_vertices.^2, 1));
parms.sbjWarp.curr_vertices = parms.sbjWarp.curr_vertices./repmat(new_radius, 3, 1) * parms.radius;

if(parms.final_unfold)
    % Unfold if necessary
    tempMesh = sbjMesh;
    tempMesh.vertices = parms.sbjWarp.curr_vertices;
    tempMesh = MARS_simpleUnfoldMesh(tempMesh, 0);
    parms.sbjWarp.curr_vertices = tempMesh.vertices;
end


if(SD_atlas.multidim)
    sbjMesh.current_data = sbjMesh.data;
else
    sbjMesh.current_data = sbjMesh.data(end,:);
end
parms.lowest_error = ComputeObjFnGivenMeshes(sbjMesh, SD_atlas, parms);

vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(sbjMesh.vertexNbors, 1)), int32(size(parms.sbjWarp.curr_vertices, 2)), ...
    int32(sbjMesh.vertexNbors), single(parms.sbjWarp.curr_vertices));
parms.harmonic_energy = sum((sqrt(vertexDistSq2Nbors(:)) - sqrt(sbjMesh.vertexDistSq2Nbors(:))).^2)./sum(sum(sbjMesh.vertexNbors~=0));
    
disp(['Orig Error: ' num2str(orig_objfn) ', Final Error: ' num2str(parms.lowest_error) ', Harm Energy: ' num2str(parms.harmonic_energy)]);


if(parms.movie_flag && ~SD_atlas.multidim)
    warped_mean = MARS_linearInterpolate(parms.sbjWarp.curr_vertices, SD_atlas.basic_atlas{level}, SD_atlas.basic_atlas{level}.mean);
    
    warped_mean2 = warped_mean; warped_mean2(1)=max_value; warped_mean2(2) = min_value;
    figure; set(gcf, 'renderer', 'zbuffer'); TrisurfMeshData(sbjMesh, warped_mean2); shading interp; colorbar; view(250, 20); 
    title('Final atlas mean after interpolating the final warp'); saveas(gcf, ['temp' num2str(parms.curr_level*4+1, '%.3d') '.tiff']);
    close all;
end











