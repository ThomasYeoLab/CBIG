function [mean_V1_perc, mean_CS_hd] = CBIG_hMRF_compute_architectonic_metrics(label, mesh_type, V1_anchors,...
     CS_anchors, num_cluster_per_hemi, debug_out_dir)
% [mean_V1_perc, mean_CS_hd] = CBIG_hMRF_compute_architectonic_metrics(label, mesh_type, V1_anchors, CS_anchors,...
% num_cluster_per_hemi, debug_out_dir)
%
% This function computes how well the current parcellation is aligned to the central sulcus and V1 boundaries.
 
% For the notations below:
% M = no of cortical vertices for both hemispheres;

% Input
%   - label: (matrix)
%     Mx1 resultant parcellation label at the current step.
%   - mesh_type: (string)
%     The string indicating the desired mesh type.
%   - V1_anchors: (character array)
%     Path to the file that stores the binary array indicating the location for V1.
%   - CS_anchors: (character array)
%     Path to the file that stores the binary array indicating the location for central sulcus.
% Output
%   - mean_V1_perc (double)
%     Percentage of the # of vertices of V1 parcels over # of total vertices within architectonic boundary of V1.
%   - mean_CS_hd (double)
%     Distance between the division between motor cortex parcels and somatosensory parcels,
%     and the architectonic boundary between area 3 and 4.
%
% Example
%   - [mean_V1_perc, mean_CS_hd] = CBIG_hMRF_compute_architectonic_metrics(label, mesh_type,...
%     V1_anchors, CS_anchors, num_cluster_per_hemi, debug_out_dir)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

[~, ~, lh_mask, rh_mask]  = CBIG_hMRF_load_mesh_mask_by_mesh_type(mesh_type);
[lh_full_label, rh_full_label] = CBIG_hMRF_get_left_right_overlapping_labels(lh_mask, rh_mask,...
    label, num_cluster_per_hemi);

folder_to_write = fullfile(debug_out_dir, ['temp_folder' num2str(rand(1))]);
[lh_full_label, rh_full_label] = CBIG_project_fsaverage2fsLR(lh_full_label, rh_full_label,...
    mesh_type, 'label', folder_to_write);

[lh_perc, rh_perc] = check_if_v1_is_vertically_striped(lh_full_label, rh_full_label, V1_anchors);
mean_V1_perc = (lh_perc + rh_perc)/2;
[hausdorff_lh, hausdorff_rh] = compute_hausdorff_dist_to_central_sulcus(lh_full_label, rh_full_label, CS_anchors);
mean_CS_hd = (hausdorff_lh + hausdorff_rh)/2;
end

function [hd_lh, hd_rh] = compute_hausdorff_dist_to_central_sulcus(lh_full_label, rh_full_label, CS_anchors)

load(CS_anchors, 'lh_CS', 'rh_CS'); 

lh_avg_mesh = CBIG_read_fslr_surface('lh', 'fs_LR_32k', 'very_inflated','medialwall.annot');
rh_avg_mesh = CBIG_read_fslr_surface('rh', 'fs_LR_32k', 'very_inflated','medialwall.annot');

[lh_full_labels,rh_full_labels] = BuildTwoVertThickBoundary(lh_avg_mesh, rh_avg_mesh,...
    lh_full_label, rh_full_label); 
lh_boundary_verts = (lh_full_labels == 0);
rh_boundary_verts = (rh_full_labels == 0);

[hd_lh, hd_rh] = compute_hausdorff_distance(lh_avg_mesh,rh_avg_mesh, lh_CS, rh_CS,...
    lh_boundary_verts,rh_boundary_verts);
end

function [lh_up_labels,rh_up_labels] = BuildTwoVertThickBoundary(lh_avg_mesh, rh_avg_mesh,...
    lh_full_label, rh_full_label)   
% this function defines the two-vertex thick boundary for a given parcellation

if(size(lh_full_label, 1) == 1)
    lh_full_label = lh_full_label';
end

if(size(rh_full_label, 1) == 1)
    rh_full_label = rh_full_label';
end

lh_neigh_labels = lh_avg_mesh.vertexNbors;
rh_neigh_labels = rh_avg_mesh.vertexNbors;
lh_neigh_labels(lh_neigh_labels ~= 0) = lh_full_label(lh_neigh_labels(lh_neigh_labels ~= 0));
rh_neigh_labels(rh_neigh_labels ~= 0) = rh_full_label(rh_neigh_labels(rh_neigh_labels ~= 0));

lh_full_labels = [lh_full_label'; lh_neigh_labels];
rh_full_labels = [rh_full_label'; rh_neigh_labels];
lh_temp = bsxfun(@minus, lh_full_labels, lh_full_labels(1,:));
rh_temp = bsxfun(@minus, rh_full_labels, rh_full_labels(1,:));
lh_temp(lh_full_labels == 0) = 0;
rh_temp(rh_full_labels == 0) = 0;
lh_up_labels = lh_full_label;
rh_up_labels = rh_full_label;
lh_up_labels(sum(lh_temp) ~= 0) = 0;
rh_up_labels(sum(rh_temp) ~= 0) = 0;
end


function [hausdorff_lh,hausdorff_rh] = compute_hausdorff_distance(lh_avg_mesh, rh_avg_mesh,...
    lh_areal_boundary_verts, rh_areal_boundary_verts, all_lh_boundary, all_rh_boundary)
% areal boundary verts corresponds to the boundary verts of a single target area based on a given parcellation
% all boundary verts corresponds to the boundary verts of all areas based on the reference segmentation
hausdorff_lh = compute_hausdorff_distance_per_hemi(lh_avg_mesh, lh_areal_boundary_verts', all_lh_boundary);
hausdorff_rh = compute_hausdorff_distance_per_hemi(rh_avg_mesh, rh_areal_boundary_verts', all_rh_boundary); 
end

function [distVec, distVecOverall] = compute_hausdorff_distance_per_hemi(hemi_avg_mesh, areal_boundary_verts,...
    all_areas_boundary_verts)
% compute Hausdorff distance per hemisphere
vertexDistSq2Nbors = MARS_computeVertexDistSq2Nbors(int32(size(hemi_avg_mesh.vertexNbors, 1)),...
    int32(size(hemi_avg_mesh.vertices, 2)), int32(hemi_avg_mesh.vertexNbors), single(hemi_avg_mesh.metricVerts));
vertexDist2Nbors = sqrt(vertexDistSq2Nbors);
[sparse_distance]=compute_sparse_distance(hemi_avg_mesh, vertexDist2Nbors);

if(~isempty(find(areal_boundary_verts == 1, 1)) && ~isempty(find(all_areas_boundary_verts == 1, 1)))
    objVec1 = (areal_boundary_verts == 1); 
    objVec2 = (all_areas_boundary_verts == 1);
    distance = compute_distance(objVec1, objVec2, sparse_distance);
    distVec(1,:) = distance;
else
    distVec(1,:) = -1;
end

distVecOverall = mean(distVec(distVec ~= -1));
end

function sparse_distance = compute_sparse_distance(avg_mesh, vertexDist2Nbors)
% compute distance between each pairs of vertices and store as sparse matrix
vertices = max(size(avg_mesh.vertexNbors));
r = reshape(repmat(13:vertices,6,1), 1, 6*(vertices-12));
sparse_distance = sparse(double(r), double(reshape(avg_mesh.vertexNbors(1:6,13:end), [1,6*(vertices-12)])),... 
double(reshape(vertexDist2Nbors(:,13:end), (size(r)))));
for i = 1:12
    sparse_distance(i,avg_mesh.vertexNbors(1:5,i)) = vertexDist2Nbors(1:5, i);
end
end

function distance = compute_distance(boundaryVec1, boundaryVec2, sparse_distance)
% compute distance between two given boundaries on the surface mesh
boundary = find(boundaryVec1);
distance_vec = zeros(size(boundaryVec1));
for idx = 1:length(boundary)
    [d, ~] = shortest_paths(sparse_distance, boundary(idx));
    distance_vec(boundary(idx)) = min(d(boundaryVec2));
end
distance = mean(nonzeros(distance_vec));
end

function [lh_perc_anchor_parcel_in_V1, rh_perc_anchor_parcel_in_V1] = check_if_v1_is_vertically_striped(lh_label,...
     rh_label, V1_anchors)

% check if parcels are vertically arranged with in area V1

load(V1_anchors,'lh_V1_vertices', 'rh_V1_vertices', 'anchors_upper_V1');

% lh
lh_anchor_parcels = lh_label(anchors_upper_V1);
lh_anchor_parcel_verts = ismember(lh_label, lh_anchor_parcels);
lh_anchor_parcel_verts = find(lh_anchor_parcel_verts == 1);
lh_perc_anchor_parcel_in_V1 = size(intersect(lh_anchor_parcel_verts, lh_V1_vertices))/size(lh_V1_vertices);
% rh
rh_anchor_parcels = rh_label(anchors_upper_V1);
rh_anchor_parcel_verts = ismember(rh_label, rh_anchor_parcels);
rh_anchor_parcel_verts = find(rh_anchor_parcel_verts == 1);
rh_perc_anchor_parcel_in_V1 = size(intersect(rh_anchor_parcel_verts, rh_V1_vertices))/size(rh_V1_vertices);
end