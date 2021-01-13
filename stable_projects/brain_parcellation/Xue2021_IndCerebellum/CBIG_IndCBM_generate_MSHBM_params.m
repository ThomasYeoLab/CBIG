function clustered = CBIG_IndCBM_generate_MSHBM_params(lh_profile, rh_profile, lh_labels, rh_labels)

% clustered = CBIG_IndCBM_generate_MSHBM_params(lh_profile, rh_profile, lh_labels, rh_labels)
%
% This function can generate parameters used in MSHBM algorithm for a given
% group-level parcellation with average profiles. 
% Input labels should match the resolution of profile. Medial wall should 
% be labeled as zero.
%
% Input:
%     - lh_profile: 
%           Nifti file of average lh profile file of your dataset. Could be
%           generated using CBIG_IndCBM_compute_profile.sh
%           Example:
%           'MSHBM_dir/profiles/avg_profile/lh_fsaverage5_roifsaverage3_avg_profile.nii.gz'
%
%     - rh_profile: 
%           Nifti file of average rh profile file of your dataset. Could be
%           generated using CBIG_IndCBM_compute_profile.sh
%           Example:
%           'MSHBM_dir/profiles/avg_profile/rh_fsaverage5_roifsaverage3_avg_profile.nii.gz'
% 
%     - lh_labels: 
%           N x 1 vector. Group-level cerebral cortical parcellation labels of lh.
% 
%     - rh_labels: 
%           N x 1 vector. Group-level cerebral cortical parcellation labels of rh.
% 
% Example:
% clustered = CBIG_IndCBM_generate_MSHBM_params('profile_dir/lh_avg_profile.nii.gz', ...
% 'profile_dir/lh_avg_profile.nii.gz', lh_labels, rh_labels)
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Check input arguments
if(size(lh_labels, 2) > 1)
    error('Input argument ''lh_labels'' should be a column vector.');
end
if(size(rh_labels, 2) > 1)
    error('Input argument ''rh_labels'' should be a column vector.');
end

lh_series = read_fmri(lh_profile);
rh_series = read_fmri(rh_profile);
series = [lh_series; rh_series];
clear lh_series rh_series
labels = [lh_labels; rh_labels];

% Check resolution
if(size(series, 1) ~= length(series))
    error('Input resolution of profile and labels do not match.');
end
D = size(series, 2) - 1;

% 0 label and zero correlation will be masked
medial_mask = labels ~= 0;
non_zero_corr_index = (sum(series, 2) ~= 0);
mask = medial_mask & non_zero_corr_index;

% Compute parameters based on labels
r = zeros(length(labels), max(labels));
rowindx = [1:1:length(labels)]';
r(sub2ind(size(r), rowindx(labels ~= 0), labels(labels ~= 0))) = 1;
r = r(mask, :);
x = series(mask, :);
x = bsxfun(@minus, x, mean(x, 2));
x = bsxfun(@times, x, 1./sqrt(sum(x.^2, 2)));
mtc = x' * r;
mtc = bsxfun(@times, mtc, 1./sqrt(sum((mtc).^2)));

clustered.d = D;
clustered.r = r;
clustered.mtc = mtc';

end

function vol = read_fmri(fmri_name)

% [fmri, vol] = read_fmri(fmri_name)
% Given the name of averaged correlation profile file (fmri_name), this
% function read in the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%       It can be .nii.gz file for correlation profile in fsaverage*
%       spaces; or .mat file for correlation profile in fs_LR_32k space
%       (there must be a variable called 'profile_mat' in the .mat file).
%
% Output:
%     - vol:
%       A num_voxels x num_timepoints matrix which is the reshaped 'vol'
%       structure for .nii.gz file or the variable 'profile_mat' for .mat file.
%

if (~isempty(strfind(fmri_name, '.nii.gz')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = fmri.vol;
    vol_size = size(vol);
    if(length(vol_size) < 4)
        vol = reshape(vol, prod(vol_size(1:3)), 1);
    else
        vol = reshape(vol, prod(vol_size(1:3)), vol_size(4));
    end
    fmri.vol = [];
else
    load(fmri_name);
    vol = profile_mat;
end

end
