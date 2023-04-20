function mismatched_parcels = CBIG_hMRF_find_parcels_mismatched_topology_indirect_nbors(lh_full_label,...
     rh_full_label, threshold, expansion_steps, mesh_type)
% mismatched_parcels = CBIG_hMRF_find_parcels_mismatched_topology_indirect_nbors(lh_full_label,...
% rh_full_label, threshold, expansion_steps, mesh_type)
% 
% This function looks for homotopic parcel pairs with different neighbors.

% For the notations below:
% N = no of vertices per hemisphere; 

% Input
%   - lh_full_label: (matrix)
%     Nx1 left hemisphere parcellation label, including medial wall vertices (marked as zeros).
%   - rh_full_label: (matrix)
%     Nx1 right hemisphere parcellation label, including medial wall vertices (marked as zeros).
%   - threshold: (double)
%     Given a threshold x, we consider two parcels to be homotopic if more than x% of their neighbors are the same.
%   - expansion_steps: (integer)
%     Say we take x expansion steps. We will take x steps' walk from the boundary of a parcel, and consider
%     all the parcels encountered during the walk as a neighbor.
%   - mesh_type: (string)
%     Now only supports 'fsaverage6'.
% Output
%   - mismatched_parcels: (matrix)
%     A list of parcel pairs with mismatched neighbors.
%
% Example
%   - mismatched_parcels = CBIG_hMRF_find_parcels_mismatched_topology_indirect_nbors(lh_full_label, rh_full_label)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~exist('threshold', 'var'))
    threshold = 0.4;
end

if(~exist('expansion_steps', 'var'))
    expansion_steps = 3;
end

if(~exist('mesh_type', 'var'))
    mesh_type = 'fsaverage6';
end

[lh_avg_mesh, rh_avg_mesh, ~, ~]  = CBIG_hMRF_load_mesh_mask_by_mesh_type(mesh_type);
if(isempty(nonzeros(intersect(lh_full_label, rh_full_label))))
    rh_full_label(rh_full_label~=0) = rh_full_label(rh_full_label~=0) - max(lh_full_label);
end

mismatched_parcels = zeros(1, max(lh_full_label));
for parcel_idx = 1:max(lh_full_label)
    lh_nearby_parcels = find_nearby_parcels(parcel_idx, lh_full_label, lh_avg_mesh, expansion_steps);
    lh_cur_nbor_set = unique(lh_nearby_parcels,'sorted');
    
    rh_nearby_parcels = find_nearby_parcels(parcel_idx, rh_full_label, rh_avg_mesh, expansion_steps);
    rh_cur_nbor_set = unique(rh_nearby_parcels,'sorted');

    % for each pair of homotopic parcels, say one has neighbors set A, another has neighbors set B
    % we consider this pair as matched if A and B are identical
    if(isequal(lh_cur_nbor_set, rh_cur_nbor_set)) 
        continue;
    % we also consider this pair as matched if A belongs to B or B belongs to A
    elseif(all(ismember(lh_cur_nbor_set, rh_cur_nbor_set)) || all(ismember(rh_cur_nbor_set, lh_cur_nbor_set))) 
        continue;
    elseif(sum(lh_cur_nbor_set == 0) ~= sum(rh_cur_nbor_set == 0))
        mismatched_parcels(parcel_idx) = 1;
    elseif(length(intersect(lh_cur_nbor_set, rh_cur_nbor_set))/length(union(lh_cur_nbor_set, rh_cur_nbor_set))...
            < threshold)
        mismatched_parcels(parcel_idx) = 1;
    else
        continue;
    end
end
mismatched_parcels = find(mismatched_parcels == 1);
end


function [nearby_parcels] = find_nearby_parcels(target_parcel_label, hemi_full_label, hemi_avg_mesh, expansion_steps)
% find nearby parcels by walking a number of 'expansion_steps' from the current target parcel boundary
full_neighborhood = hemi_avg_mesh.vertexNbors;
target_parcel_verts = find(target_parcel_label == hemi_full_label);
random_vert_from_cur_parcel = target_parcel_verts(1);
if(length(target_parcel_verts) < 7)
    boundary_verts = target_parcel_verts;
else
    boundary_verts=find_parcel_boundary_including_medial_wall(hemi_avg_mesh, hemi_full_label, target_parcel_label);
end

cur_layer_verts_idx = boundary_verts;
nearby_parcels = [];
is_adjacent_to_medial = 0;
for cur_layer_idx = 1:expansion_steps
    cur_nborhood = full_neighborhood(:, cur_layer_verts_idx);
    cur_nborhood(cur_nborhood == 0) = random_vert_from_cur_parcel;
    cur_nborhood_labels = hemi_full_label(cur_nborhood);
    outer_layer_verts_idx = unique(cur_nborhood(cur_nborhood_labels ~= target_parcel_label));

    if(sum(cur_nborhood_labels(:) == 0) && cur_layer_idx == 1)
        is_adjacent_to_medial = 1;         
    end       

    nearby_parcels = [nearby_parcels unique(hemi_full_label(outer_layer_verts_idx))];
    hemi_full_label(outer_layer_verts_idx) = target_parcel_label;
    cur_layer_verts_idx = outer_layer_verts_idx;
end

nearby_parcels = unique(nearby_parcels);
if(is_adjacent_to_medial)
    nearby_parcels = [nearby_parcels 0];
end
end


function [boundary_verts] = find_parcel_boundary_including_medial_wall(hemi_avg_mesh, full_label, target_parcel)
% this function detects the boundary vertices for the given parcel based on the given hemisphere
target_parcel_verts = find(full_label == target_parcel);

neighborhood = hemi_avg_mesh.vertexNbors;
random_vert_from_target_parcel = target_parcel_verts(1);
target_neighborhood = neighborhood(:, target_parcel_verts);
target_neighborhood(target_neighborhood == 0) = random_vert_from_target_parcel;

neighborhood_label = full_label(target_neighborhood);
neighborhood_label_is_target = neighborhood_label == target_parcel;
boundary_verts = target_parcel_verts(sum(neighborhood_label_is_target,1) ~= 6);
end

