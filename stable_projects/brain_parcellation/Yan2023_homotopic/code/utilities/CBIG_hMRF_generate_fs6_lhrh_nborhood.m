function lh_rh_Neighborhood = CBIG_hMRF_generate_fs6_lhrh_nborhood(outdir)
% CBIG_hMRF_generate_fs6_lhrh_nborhood(outdir)
% 
% The left and right fsaverage6 hemispheres were first registered using the Freesurfer's 'mris_left_right_register'
% function. After registration, the closest right hemisphere vertex was found for each left hemisphere vertex.
% Similarly, the closest left hemisphere vertex was found for each right hemisphere vertex. Here, the closest vertex
% was defined based on geodesic distance on the spherical surface mesh. A pair of left and right hemisphere vertices
% was considered homotopic if they were each other's closest neighbor. 
%
% Optional input
%   - output_dir: (string) 
%     If given, it is the ABSOLUTE path to the directory to which the output results will be saved.
%
% Output
%   - lh_rh_Neighborhood: (matrix)
%     This is the output cross-hemispheric neighborhood where by a left hemisphere vertex X and a right
%     hemisphere vertex Y are considered neighbors iff they are closest to one another (as explained above).
%
% Example
%   - CBIG_hMRF_generate_fs6_lhrh_nborhood(your_outdir)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
hMRF_code_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Yan2023_homotopic', 'code');
fs6_surf_dir = fullfile(hMRF_code_dir, 'utilities', 'input'); % the folder contains freesurfer fsaverage6 info

% match by min distance
lh_mesh6_sym = ReadNCAvgMesh_additional_input_var('lh', 'fsaverage6', 'sphere.left_right', 'cortex', fs6_surf_dir);
rh_mesh6_sym = ReadNCAvgMesh_additional_input_var('rh', 'fsaverage6', 'sphere.left_right', 'cortex', fs6_surf_dir);
lh_verts = lh_mesh6_sym.vertices;
rh_verts = rh_mesh6_sym.vertices;

% find pairwise euclidean distance
dist_mat = pdist2(lh_verts', rh_verts', 'euclidean');

% find the closest rh vert for each lh vert
rh2lh_assigned_matched_vert = zeros(size(dist_mat,1),1);
for row_idx = 1:size(dist_mat,1)
    [~, rh2lh_assigned_matched_vert(row_idx)] = min(dist_mat(row_idx,:));
end

% find the closest lh vert for each rh vert
lh2rh_assigned_matched_vert = zeros(size(dist_mat,1),1);
for col_idx = 1:size(dist_mat,1)
    [~, lh2rh_assigned_matched_vert(col_idx)] = min(dist_mat(:,col_idx));
end

% find lh_one2rh_one pairs
lh_one2rh_one = zeros(size(dist_mat,1),1);
lh_many_to_rh_one_counter = 0;
lh_one_to_rh_one = 0;
lh_one_to_rh_many_counter = 0;
lh_not_matched = 0;

for lh_vert_idx = 1:length(rh2lh_assigned_matched_vert)
    min_dist_rh_vert = rh2lh_assigned_matched_vert(lh_vert_idx);
    if(sum(rh2lh_assigned_matched_vert == min_dist_rh_vert) == 1) 
        % means that cur rh vert is matched to 1 lh vert
        if(lh2rh_assigned_matched_vert(min_dist_rh_vert) == lh_vert_idx)  
            % means that cur lh vert is matched to the given rh vert
            if(sum(lh2rh_assigned_matched_vert == lh_vert_idx) == 1) 
                % means that cur lh vert is matched to 1 rh vert
                lh_one2rh_one(lh_vert_idx) = min_dist_rh_vert;
                lh_one_to_rh_one = lh_one_to_rh_one + 1;
            else % means that cur lh vert is matched to many rh vert
                lh_one_to_rh_many_counter = lh_one_to_rh_many_counter + 1;
            end
        else % means that cur lh vert is matched to a rh vert, which does not match back
            lh_not_matched = lh_not_matched + 1;
        end
    else % means that cur rh vert is matched to more than one lh vert
        lh_many_to_rh_one_counter = lh_many_to_rh_one_counter + 1;
    end
end

% validation on lh
lh_sum = lh_many_to_rh_one_counter + lh_one_to_rh_one + lh_one_to_rh_many_counter + lh_not_matched;

% make the one-2-one vertices pairs into a NxN matrix (N = no. vertices per hemisphere)
lh_rh_Neighborhood = zeros(length(lh_one2rh_one), length(lh_one2rh_one));
for lh_vert_idx = 1:length(lh_one2rh_one)
    if(lh_one2rh_one(lh_vert_idx)) % this lh vert is matched to some vert on rh
        lh_rh_Neighborhood(lh_vert_idx, lh_one2rh_one(lh_vert_idx)) = 1;
    end
end

if(exist('outdir', 'var'))
    save(fullfile(outdir, 'fs_nborhood_derived_from_fs6.mat'), 'lh_rh_Neighborhood', '-v7.3');
end

end


function avg_mesh = ReadNCAvgMesh_additional_input_var(hemi, mesh_name, surf_type, label, SUBJECTS_DIR)

% Read the average mesh (from FREESURFER) of normal control people. This is modified based on CBIG_ReadNCAvgMesh
% since here we need to specify SUBJECTS_DIR.
% 
%   avg_mesh = ReadNCAvgMesh_additional_input_var(hemi, mesh_name, surf_type, label, SUBJECTS_DIR)
%   Input:
%       hemi        : 'lh' or 'rh'
%       mesh_name   : 'fsaverage6', 'fsaverage5', 'fsaverage4'
%       surf_type   : 'inflated', 'white', 'sphere'
%       label       : 'cortex', 'Yeo2011_7Networks_N1000.annot'
%   Output: 
%       avg_mesh    : average mesh
% 
%   PS: For more labels, check $FREESURFER_HOME/subjects/fsaverage/label
%
% Written by Xiaoxuan Yan CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(nargin < 3)
    surf_type = 'white'; 
end

if(~exist('SUBJECTS_DIR','var'))
    SUBJECTS_DIR = fullfile(getenv('FREESURFER_HOME'), 'subjects');
end
parms.SUBJECTS_DIR = SUBJECTS_DIR;
parms.hemi = hemi;
parms.SUBJECT = mesh_name;
parms.read_surface = @MARS2_readSbjMesh;
parms.radius = 100;
parms.unfoldBool = 0;
parms.flipFacesBool = 1;
parms.surf_filename = 'sphere';
parms.metric_surf_filename = surf_type;
parms.data_filename_cell = {'sulc', 'curv'};

if(nargin == 5)
    if(strfind(label, 'annot'))
        parms.annot_filename = label;
    else
        parms.label_filename_cell = {label};
    end
end

avg_mesh = MARS2_readSbjMesh(parms);
avg_mesh.vertices = avg_mesh.metricVerts;

% Recompute vertexDistSq2Nbors and faceAreas
avg_mesh.vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(avg_mesh.vertexNbors, 1)),...
 int32(size(avg_mesh.vertices, 2)), int32(avg_mesh.vertexNbors), single(avg_mesh.vertices));
avg_mesh.faceAreas = MARS_computeMeshFaceAreas(int32(size(avg_mesh.faces, 2)), int32(avg_mesh.faces),...
 single(avg_mesh.vertices));
end