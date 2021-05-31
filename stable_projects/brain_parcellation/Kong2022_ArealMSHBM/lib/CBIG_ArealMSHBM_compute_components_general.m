function parcel_components = CBIG_ArealMSHBM_compute_components_general(lh_labels,rh_labels,mesh,num_parcel)

% parcel_components = CBIG_ArealMSHBM_compute_components_general(lh_labels,rh_labels,mesh,num_parcel)
%
% This script will compute the number of components for each parcel of a given parcellation <lh_labels>
% <rh_labels>. If there is any empty parcel, the number of components for empty parcels will be denoted as
% NaN.
%
% Input:
%   - lh_labels, rh_labels: (Nx1 vector)
%     The given parcellation labels. Medial wall should be denoted as 0.
%
%   - mesh: (string)
%     The mesh name of the parcellation surface space. For example, 'fsaverage6', 'fs_LR_32k'.
%
%   - num_parcel: (scalar)
%     The total number of parcels (excluding medial wall). For example, 400.
%
% Output:
%   - parcel_componets: (1xnum_parcel vector)
%     The number of components for each parcel of a given parcellation.
%
% Example:
%   parcel_components = CBIG_ArealMSHBM_compute_components_general(lh_labels,rh_labels,'fs_LR_32k',400);
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md



if(strcmp(mesh, 'fs_LR_32k'))
    lh_avg_mesh=CBIG_read_fslr_surface('lh','fs_LR_32k','inflated');
    rh_avg_mesh=CBIG_read_fslr_surface('rh','fs_LR_32k','inflated');
else
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');
end

% Find empty parcels
for i = 1:num_parcel/2
    if(sum(lh_labels == i) == 0)
        lh_labels_mask(i) = 1;
    else
        lh_labels_mask(i) = 0;
    end
end

for i = (num_parcel/2 + 1):num_parcel
    if(sum(rh_labels == i) == 0)
        rh_labels_mask(i-num_parcel/2) = 1;
    else
        rh_labels_mask(i-num_parcel/2) = 0;
    end
end

% Generate sizes of components
[~, lh_sizes, ~] = CBIG_ComputeConnectedComponentsFromSurface(lh_avg_mesh, lh_labels);
[~, rh_sizes, ~] = CBIG_ComputeConnectedComponentsFromSurface(rh_avg_mesh, rh_labels);

% For non-empty parcels, find the number of components
labels_mask = [lh_labels_mask rh_labels_mask];
parcel_components(labels_mask==0) = [lh_sizes(2:end) rh_sizes(2:end)];
parcel_components = cellfun('length',parcel_components);

% For empty parcels, set the number of componnets to be NaN
parcel_components(labels_mask==1) = 0;
parcel_components(parcel_components == 0) = NaN;
