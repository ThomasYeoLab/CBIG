function hemi_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(hemi_full_corr,...
    target_parcel_list, hemi_full_label, hemi_avg_mesh, hemi_mask)
% hemi_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(hemi_full_corr,...
% target_parcel_list, hemi_full_label, hemi_avg_mesh, hemi_mask)
%
% This function merges the parcels indicated in the list `target_parcel_list` given the parcellation label.
% Singular parcel means that this parcel does not have a homotopic counterpart on the opposite hemisphere.

% For the notations below:
% N = no of vertices per hemisphere; 
% M = no of cortical vertices for both hemispheres;
% k = total no of parcels 

% Input
%   - hemi_full_corr: (struct)
%     The NxN partical premultiplied matrix for a single hemisphere.
%   - target_parcel_list: (matrix)
%     The list of parcels to be merged.
%   - hemi_full_label: (matrix)
%     Resultant Nx1 parcellation label at the current step for a single hemisphere.
%   - hemi_avg_mesh: (struct)
%     The struct containing mesh information for a single hemisphere.
%   - hemi_mask: (matrix)
%     A Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on a single hemisphere.
%
% Output
%   - hemi_full_label: (matrix)
%     Resultant Nx1 parcellation label after merging the given list of parcels.
%
% Example
%   - hemi_full_label = CBIG_hMRF_merge_singular_parcels_on_one_hemi(hemi_full_corr,...
%     target_parcel_list, hemi_full_label, hemi_avg_mesh, hemi_mask)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

full_neighborhood = hemi_avg_mesh.vertexNbors;

for target_parcel_idx = 1:length(target_parcel_list)
    target_parcel_label = target_parcel_list(target_parcel_idx);
    target_parcel_verts = find(target_parcel_label == hemi_full_label);
    random_vert_from_cur_parcel = target_parcel_verts(1);
    if(length(target_parcel_verts) < 7)
        boundary_verts = target_parcel_verts;
    else
        [boundary_verts] = find_parcel_boundary_non_medial_wall(hemi_avg_mesh, hemi_full_label,...
            target_parcel_label, hemi_mask);
    end
    layer_to_reassign_verts_idx = boundary_verts;
    
    while(true)
        cur_nborhood = full_neighborhood(:, layer_to_reassign_verts_idx);
        cur_nborhood(cur_nborhood == 0) = random_vert_from_cur_parcel;
        surrounding_nborhood_labels = find_labels_of_surrounding_nbors_for_cur_layer(cur_nborhood,...
        hemi_full_label,target_parcel_label);
        
        if(size(surrounding_nborhood_labels,1) ~= 6)
            surrounding_nborhood_labels = surrounding_nborhood_labels';
        end
        
        % Reassigning labels according to neighboring vertices...
        for vertex_order = 1:length(layer_to_reassign_verts_idx)
            vertex_idx = layer_to_reassign_verts_idx(vertex_order);
            vertex_nborhood_labels = surrounding_nborhood_labels(:,vertex_order);
            hemi_full_label = reassign_layer_vertex(vertex_nborhood_labels, hemi_full_label,...
                hemi_full_corr, vertex_idx);
        end
        
        cur_nborhood_labels = hemi_full_label(cur_nborhood);
        inner_layer_verts_idx = unique(cur_nborhood(cur_nborhood_labels == target_parcel_label));
        if(isempty(inner_layer_verts_idx))
            break;
        end
        layer_to_reassign_verts_idx = inner_layer_verts_idx;
    end
end
end

function [boundary_verts] = find_parcel_boundary_non_medial_wall(hemi_avg_mesh, full_label,...
    target_parcel, hemi_mask)
% this function detects the boundary vertices for the given parcel based on the given hemisphere

target_parcel_verts = find(full_label == target_parcel);
neighborhood = hemi_avg_mesh.vertexNbors;
random_vert_from_target_parcel = target_parcel_verts(1);
target_neighborhood = neighborhood(:, target_parcel_verts);
target_neighborhood(target_neighborhood == 0) = random_vert_from_target_parcel;

neighborhood_label = full_label(target_neighborhood);
neighborhood_label_is_target = neighborhood_label == target_parcel;
boundary_verts = target_parcel_verts(sum(neighborhood_label_is_target,1) ~= 6);
medial_wall_verts = find(hemi_mask == 0);
boundary_verts = setdiff(boundary_verts, medial_wall_verts); % need to exclude medial wall boundary  
end

function cur_nborhood_labels = find_labels_of_surrounding_nbors_for_cur_layer(cur_nborhood,...
    hemi_full_label,target_parcel_label)
% find the labels of the surrounder neighbor vertices of the current layer of vertices
cur_nborhood_labels = hemi_full_label(cur_nborhood);
cur_nborhood_labels(cur_nborhood_labels == target_parcel_label) = 0;
end

function hemi_full_label = reassign_layer_vertex(vertex_nborhood_labels, hemi_full_label,...
    hemi_full_corr, vertex_idx)
% this function reassigns a vertex from the current layer with the proper label

candidate_labels = nonzeros(unique(vertex_nborhood_labels));
if(length(candidate_labels) == 1) % there is only 1 candidate label to choose from
    hemi_full_label(vertex_idx) = candidate_labels;
elseif(length(candidate_labels) > 1) % we have competing labels to choose from
    label_to_assign = choose_which_label_to_assign(vertex_idx, candidate_labels, hemi_full_label, hemi_full_corr);
    hemi_full_label(vertex_idx) = label_to_assign; 
end
end

function label_to_assign = choose_which_label_to_assign(vertex_idx, candidate_labels,...
    hemi_full_label, hemi_full_corr)
% this function determines which label the current vertex should be assigned to
% based on maximum correlation between the current vertex and the parcel that the label corresponds to

max_corr = -1e5;
for candidate_label_idx = 1:length(candidate_labels)
    cur_parcel_label = candidate_labels(candidate_label_idx);
    cur_parcel_verts = hemi_full_label == cur_parcel_label;
    curr_corr = mean(hemi_full_corr(vertex_idx, cur_parcel_verts));
    if(mean(curr_corr(:)) > max_corr)
        max_corr = curr_corr;
        label_to_assign =  cur_parcel_label;
    else
        continue;
    end
end
end