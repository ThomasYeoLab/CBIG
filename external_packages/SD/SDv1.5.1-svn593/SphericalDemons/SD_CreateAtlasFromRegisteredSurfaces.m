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
function SD_atlas = SD_CreateAtlasFromRegisteredSurfaces(SD_atlas)

% Given a bunch of registered brains, create the brain atlas.
%
%
% SD_atlas.zrotate = pi/90;
% SD_atlas.yrotate = pi/90;
% SD_atlas.radius = 100;
% SD_atlas.basic_atlas = [];
% 
% % These parameters are mostly used for reading subjects and meshes
% 
% % parameters to be filled
% SD_atlas.parms.hemi = 'lh';
% SD_atlas.parms.SUBJECTS_DIR = '/afs/csail/group/vision/cortex/ythomas/ManualLabeled/';
% SD_atlas.parms.uniform_mesh_dir = '/afs/csail/group/vision/cortex/ythomas/work/MARS/';
% %SD_atlas.parms.warp_filename = '';
% SD_atlas.parms.subject_cell = '';
% %SD_atlas.parms.annot_filename = '';
% 
% % parameters likely to stay the same
% SD_atlas.parms.surf_filename = 'sphere';
% SD_atlas.parms.data_filename_cell = {'inflated.H','sulc','curv'};
% SD_atlas.parms.WORK_DIR = 'SD';
% SD_atlas.parms.uniform_meshes = {'ic4.tri', 'ic5.tri', 'ic6.tri', 'ic7.tri'};


if(~SD_atlas.multidim && length(SD_atlas.parms.data_filename_cell) ~= length(SD_atlas.parms.uniform_meshes))
   error('Num data not equal to number of meshes!!'); 
end

num_meshes = length(SD_atlas.parms.uniform_meshes);
num_subjects = length(SD_atlas.parms.subject_cell);
SD_atlas.basic_atlas = cell(num_meshes, 1);

temp_parms = SD_atlas.parms;
SD_sbjMeshCell = cell(num_subjects, 1);
for j = 1:num_subjects
    
    disp(['Reading subject ' num2str(j)]);
    temp_parms.SUBJECT = SD_atlas.parms.subject_cell{j};
    if(~isfield(temp_parms, 'read_surface'))
        SD_sbjMesh = MARS2_readSbjMesh(temp_parms);
    else
        SD_sbjMesh = feval(temp_parms.read_surface, temp_parms);
    end
        
    if(isfield(temp_parms, 'warp_filename'))
        disp(['Reading warp ' num2str(temp_parms.warp_filename)]);
        sbjWarp_struct = load([temp_parms.SUBJECTS_DIR '/' temp_parms.SUBJECT '/' temp_parms.WORK_DIR '/' temp_parms.warp_filename]);
        SD_sbjMesh.vertices = sbjWarp_struct.sbjWarp.curr_vertices;
    end
    SD_sbjMeshCell{j} = SD_sbjMesh;
end


% if(SD_atlas.parms.normalizeBool == 2)
%    normalization = zeros(size(SD_sbjMeshCell{1}.data, 1), num_subjects);   
%    for j = 1:num_subjects
%        normalization(:, j) = SD_sbjMeshCell{j}.normalization;
%    end
%    normalization = sum(normalization, 2)./sum(normalization~=0, 2);
% end

max_mesh = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic7.tri');
max_mesh.vertices = MARS_yrotate(MARS_zrotate(max_mesh.vertices, SD_atlas.zrotate), SD_atlas.yrotate);

var_mesh = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, 'ic4.tri');
var_mesh.vertices = MARS_yrotate(MARS_zrotate(var_mesh.vertices, SD_atlas.zrotate), SD_atlas.yrotate);

spatial2conditional = MARS_findNearestVertex(max_mesh.vertices, var_mesh);
if(SD_atlas.multidim)
    data_matrix = zeros(num_subjects, size(SD_sbjMeshCell{1}.data, 1), size(max_mesh.vertices, 2));
    for j = 1:num_subjects
        disp(num2str(j))
        data_matrix(j, :, :) = MARS_linearInterpolate(max_mesh.vertices, SD_sbjMeshCell{j}, SD_sbjMeshCell{j}.data);
    end
    
    structure_index = zeros(size(SD_sbjMeshCell{1}.data, 1), num_subjects);
    BIG_NO = -1e6;
    for j = 1:size(SD_sbjMeshCell{1}.data, 1)
        structure_index(j, :) = (squeeze(data_matrix(:, j, 1)) > BIG_NO*0.1);
    end
    structure_index = logical(structure_index);
else
    data_matrix = zeros(num_subjects, size(max_mesh.vertices, 2));    
end

for i = 1:num_meshes
    disp(['Computing basic atlas ' num2str(i)]);
    basic_atlas = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, SD_atlas.parms.uniform_meshes{i});
    basic_atlas.vertices = MARS_yrotate(MARS_zrotate(basic_atlas.vertices, SD_atlas.zrotate), SD_atlas.yrotate);

    if(~SD_atlas.multidim)
        for j = 1:num_subjects
            disp(num2str(j))
            data_matrix(j, :) = MARS_linearInterpolate(max_mesh.vertices, SD_sbjMeshCell{j}, SD_sbjMeshCell{j}.data(i, :));
        end
    end
    
    % compute mean
    if(SD_atlas.multidim)
        dense_mean = zeros(size(data_matrix,2), size(data_matrix, 3));
        for j = 1:size(SD_sbjMeshCell{1}.data, 1)
            dense_mean(j, :) = single(mean(squeeze(data_matrix(structure_index(j, :), j, :)), 1));
        end    
    else
        dense_mean = squeeze(single(mean(data_matrix, 1))); % D x N
    end
    basic_atlas.mean = dense_mean(:, 1:size(basic_atlas.vertices, 2));
    
    
    % compute variance
    if(SD_atlas.multidim)
        dense_variance = zeros(size(data_matrix,2), size(data_matrix, 3));
        for j = 1:size(SD_sbjMeshCell{1}.data, 1)
            dense_variance(j, :) = single(var(squeeze(data_matrix(structure_index(j, :), j, :)), 0, 1));
        end
        sparse_variance = zeros(size(SD_sbjMeshCell{1}.data, 1), size(var_mesh.vertices, 2));
        count = zeros(size(SD_sbjMeshCell{1}.data, 1), size(var_mesh.vertices, 2));
    else
        dense_variance = squeeze(single(var(data_matrix, 0, 1)));
        sparse_variance = zeros(1, size(var_mesh.vertices, 2));
        count = zeros(1, size(var_mesh.vertices, 2));
    end
    
    if(SD_atlas.accumulate_variance)
        for j = 1:size(max_mesh.vertices, 2)
            sparse_variance(:, spatial2conditional(j)) = sparse_variance(:, spatial2conditional(j)) + dense_variance(:, j);
            count(:, spatial2conditional(j)) = count(:, spatial2conditional(j)) + 1;
        end

        sparse_variance = sparse_variance./count;
        if(SD_atlas.relative_variance)
            for j = 1:size(sparse_variance, 1)
                sparse_variance(j, :) = max(sparse_variance(j, :), 0.1*mean(sparse_variance(j, :))); %threshold
            end
        else
            sparse_variance = max(sparse_variance, SD_atlas.relative_variance_val);
        end

        if(size(sparse_variance, 2) >= size(basic_atlas.vertices, 2))
            basic_atlas.variance = sparse_variance(:, 1:size(basic_atlas.vertices, 2));
        else
            basic_atlas.variance = MARS_linearInterpolate(basic_atlas.vertices, var_mesh, sparse_variance);
        end
    else
        sparse_variance = dense_variance(:, 1:size(basic_atlas.vertices, 2));
        if(SD_atlas.relative_variance)
            for j = 1:size(sparse_variance, 1)
                sparse_variance(j, :) = max(sparse_variance(j, :), 0.1*mean(sparse_variance(j, :))); %threshold
            end
        else
            sparse_variance = max(sparse_variance, SD_atlas.relative_variance_val);
        end
        basic_atlas.variance = sparse_variance;
    end
        
%     if(SD_atlas.multidim && SD_atlas.parms.normalizeBool == 2)
%         for j = 1:size(SD_sbjMeshCell{1}.data, 1)
%             basic_atlas.variance(j, :) = basic_atlas.variance(j, :)*normalization(j);
%         end
%     end
    
    if(SD_atlas.multidim && SD_atlas.parms.scale~=0 && SD_atlas.parms.speed~=0)
        for j = 1:size(SD_sbjMeshCell{1}.data, 1)
            max_value = max(basic_atlas.mean(j, :));
            min_value = min(basic_atlas.mean(j, :));
            z = (basic_atlas.mean(j, :) - min_value)./(max_value - min_value) - 0.5;
            basic_atlas.variance(j, :) = basic_atlas.variance(j, :) .* (SD_atlas.parms.scale * 1./(1 + exp(SD_atlas.parms.speed*z)) + 1);
        end
    end
    SD_atlas.basic_atlas{i} = basic_atlas;
end


% for i = 1:num_meshes
%     SD_atlas.basic_atlas{i}.variance(SD_atlas.basic_atlas{i}.mean >= 0) = 1;
%     SD_atlas.basic_atlas{i}.variance(SD_atlas.basic_atlas{i}.mean < 0) = 10;
% end
    
    
% if(temp_parms.normalizeBool == 2)
%     for i = 1:num_meshes
%         SD_atlas.basic_atlas{i}.mean = SD_atlas.basic_atlas{i}.mean - repmat(max(SD_atlas.basic_atlas{i}.mean, [], 2), 1, size(SD_atlas.basic_atlas{i}.mean, 2)) + 1;
%     end
% end

















        
        
        
        
        
        