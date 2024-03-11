function TC = CBIG_pFIC_generate_fslr_TC_Desikan(input_path)

% TC = CBIG_pFIC_generate_fslr_TC_Desikan(input_path)
% This function generates parcellated time course (Desikan parcellation, to
% be more specific) from vertex-level time course in FSLR space. ROI
% definitions can be found in Desikan_72.txt
% Input:
%   - input_path: absolute path to the .dtseries.nii file (contains a 96854-by-T matrix, T is the number of frames)
% Output:
%   - TC: parcellated time course (72-by-T, T is the number of framas)
% Example:
% TC_100106 = CBIG_pFIC_generate_fslr_TC_Desikan('HCP/100106/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
%
% Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Load pacellation labels
% This function can be found here in our GitHub repo:
% ThomasYeoLab/CBIG/utilities/matlab/fslr_matlab/CBIG_read_fslr_surface.m
lh_desikan_fslr_32k = CBIG_read_fslr_surface('lh', 'fs_LR_32k', 'very_inflated', 'aparc.annot');
lh_desikan_fslr_32k_label = lh_desikan_fslr_32k.MARS_label;
rh_desikan_fslr_32k = CBIG_read_fslr_surface('rh', 'fs_LR_32k', 'very_inflated', 'aparc.annot');
rh_desikan_fslr_32k_label = rh_desikan_fslr_32k.MARS_label;

TC_lh = zeros(36, 1200);
TC_rh = zeros(36, 1200);

% Read cifti surface data
surface = ft_read_cifti(input_path);
surface = surface.dtseries;
surface_lh = surface(1:32492, :);
surface_rh = surface(32493:64984, :);  

% Parcellate
for i = 1:36
    TC_lh(i, :) = nanmean(surface_lh(lh_desikan_fslr_32k_label == i, :));
    TC_rh(i, :) = nanmean(surface_rh(rh_desikan_fslr_32k_label == i, :)); 
end
TC = [TC_lh; TC_rh];

end
