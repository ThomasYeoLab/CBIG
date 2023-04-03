function [lh_full_label, rh_full_label] = CBIG_hMRF_get_left_right_overlapping_labels(lh_mask, rh_mask,...
     label, num_cluster_per_hemi)
% [lh_full_label, rh_full_label] = CBIG_hMRF_get_left_right_overlapping_labels(lh_mask, rh_mask,...
% label, num_cluster_per_hemi)
%
% This function reformats the input label vector into 2 vectors for different hemispheres.
% Basically, it make left and right labels share overlapping indices, e.g., If the original parcellation 
% has labels 1~200 on the left and 201~400 on the right, this function reformats the labels into
% left 1~200 and right also 1~200.
% The resultant vectors also takes into account medial wall vertices as indicated by zeros. 

% For the notations below:
% N = no of vertices per hemisphere; 
% M = no of cortical vertices for both hemispheres;

% Input
%   - lh_mask: (matrix)
%     A Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the left hemisphere.
%   - rh_mask: (matrix)
%     A Nx1 binary array indicating whether a vertex belongs to cortex or the medial wall on the right hemisphere.
%   - label: (matrix)
%     Resultant Mx1 parcellation label at the current step.
%   - num_cluster_per_hemi: (integer)
%     The number of clusters per hemisphere.
% Output
%   - lh_full_label: (matrix)
%     Resultant Nx1 left hemisphere parcellation label, including medial wall vertices (marked as zeros).
%   - rh_full_label: (matrix)
%     Resultant Nx1 right hemisphere parcellation label, including medial wall vertices (marked as zeros).
%
% Example
% [lh_full_label, rh_full_label] = CBIG_hMRF_get_left_right_overlapping_labels(lh_mask, rh_mask,...
% label, num_cluster_per_hemi)
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% get full labels on lh and rh   
lh_label = label(1: sum(lh_mask));
rh_label = label(sum(lh_mask)+1: end);

lh_full_label(lh_mask) = lh_label;
rh_full_label(rh_mask) = rh_label;

%% make lh rh parcels overlap
if(isempty(nonzeros(intersect(lh_full_label, rh_full_label))))
    rh_full_label(rh_full_label~=0) = rh_full_label(rh_full_label~=0) - cast(num_cluster_per_hemi,...
        class(rh_full_label));
end
end