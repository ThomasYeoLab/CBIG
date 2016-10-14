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
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================
function SD_atlas = SD_NormalizeAtlasWarps(SD_atlas)

% Given a bunch of registered brains, normalize the corresponding warps so
% that the atlas coordinate system is unbiased, i.e, the average
% deformation from the atlas coordinates is zero.
%
%
% SD_atlas.zrotate = pi/90;
% SD_atlas.yrotate = pi/90;
% SD_atlas.radius = 100;
% SD_atlas.basic_atlas = [];
% 
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
temp_parms = SD_atlas.parms;

if(~isfield(temp_parms, 'warp_filename'))
    return;
end

num_subjects = length(SD_atlas.parms.subject_cell);


SD_sbjMeshCell = cell(num_subjects, 1);

basic_atlas = MARS_readUniformMesh(SD_atlas.parms.uniform_mesh_dir, SD_atlas.parms.uniform_meshes{end});
basic_atlas.vertices = MARS_yrotate(MARS_zrotate(basic_atlas.vertices, SD_atlas.zrotate), SD_atlas.yrotate);


for j = 1:num_subjects
    
    disp(['Reading subject ' num2str(j)]);
    temp_parms.SUBJECT = SD_atlas.parms.subject_cell{j};
    if(~isfield(temp_parms, 'read_surface'))
        SD_sbjMesh = MARS2_readSbjMesh(temp_parms);
    else
        SD_sbjMesh = feval(temp_parms.read_surface, temp_parms);
    end
    
    disp(['Reading warp ' num2str(temp_parms.warp_filename)]);
    sbjWarp_struct = load([temp_parms.SUBJECTS_DIR '/' temp_parms.SUBJECT '/' temp_parms.WORK_DIR '/' temp_parms.warp_filename]);
    SD_sbjMesh.orig_vertices = SD_sbjMesh.vertices;
    SD_sbjMesh.vertices = sbjWarp_struct.sbjWarp.curr_vertices;
    SD_sbjMeshCell{j} = SD_sbjMesh;
end

SD_InvsbjMeshCell = cell(num_subjects, 1);

for j = 1:num_subjects
    SD_InvsbjMeshCell{j}.vertices = basic_atlas.vertices;
    for ii = 1:3
        SD_InvsbjMeshCell{j}.vertices(ii,:) = MARS_linearInterpolate(basic_atlas.vertices, SD_sbjMeshCell{j}, SD_sbjMeshCell{j}.orig_vertices(ii,:));
    end
    SD_InvsbjMeshCell{j}.vertices = SD_InvsbjMeshCell{j}.vertices.*SD_atlas.radius./repmat(sqrt(sum(SD_InvsbjMeshCell{j}.vertices.^2,1)), [3 1]);   
end


for tmp_cnt = 1:10
    
    SD_InvsbjMeshAvg.vertices = zeros(size(SD_InvsbjMeshCell{1}.vertices));
    for j = 1:num_subjects
    
        SD_InvsbjMeshAvg.vertices = SD_InvsbjMeshAvg.vertices + (SD_InvsbjMeshCell{j}.vertices - basic_atlas.vertices);
    end
    
    SD_InvsbjMeshAvg.vertices = SD_InvsbjMeshAvg.vertices/num_subjects;
    
    display(['Average correction: ' num2str(mean(sqrt(sum(SD_InvsbjMeshAvg.vertices.^2, 1))))]);
    for j = 1:num_subjects
        SD_InvsbjMeshCell{j}.vertices = SD_InvsbjMeshCell{j}.vertices - SD_InvsbjMeshAvg.vertices;
        SD_InvsbjMeshCell{j}.vertices = SD_InvsbjMeshCell{j}.vertices.*SD_atlas.radius./repmat(sqrt(sum(SD_InvsbjMeshCell{j}.vertices.^2,1)), [3 1]);
    end
        
end


for j = 1:num_subjects
    
    tmp_mesh = basic_atlas;
    tmp_mesh.vertices = SD_InvsbjMeshCell{j}.vertices;
    temp_parms.SUBJECT = SD_atlas.parms.subject_cell{j};
    
    sbjWarp_struct = load([temp_parms.SUBJECTS_DIR '/' temp_parms.SUBJECT '/' temp_parms.WORK_DIR '/' temp_parms.warp_filename]);
    
    for ii = 1:3
        sbjWarp_struct.sbjWarp.curr_vertices(ii,:) = MARS_linearInterpolate(SD_sbjMeshCell{j}.orig_vertices, tmp_mesh, basic_atlas.vertices(ii,:));
    end
    
    sbjWarp.curr_vertices = sbjWarp_struct.sbjWarp.curr_vertices.*SD_atlas.radius./repmat(sqrt(sum(sbjWarp_struct.sbjWarp.curr_vertices.^2,1)), [3 1]);
    save([temp_parms.SUBJECTS_DIR '/' temp_parms.SUBJECT '/' temp_parms.WORK_DIR '/' temp_parms.warp_filename], 'sbjWarp');
    
end




