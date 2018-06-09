function [fmri, vol, vol_size] = CBIG_MSHBM_read_fmri(fmri_name)

% [fmri, vol, vol_size] = CBIG_MSHBM_read_fmri(fmri_name)
% Given the name of functional MRI file (fmri_name), this function read in
% the fmri structure and the content of signals (vol) or functional 
% connectivity matrix.
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%
% Output:
%     - fmri:
%       The structure read in by MRIread() or ft_read_cifti(), or a simple
%       .mat file contain the functional connectivity profile matrix
%       profile_mat. To save the memory, fmri.vol (for NIFTI) or 
%       fmri.dtseries (for CIFTI) is set to be empty after it is transfered 
%       to "vol".
%
%     - vol:
%       A num_voxels x num_timepoints matrix which is the content of
%       fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) after reshape or
%       A num_voxels x num_ROIs matrix which is the functional
%       connectivity.
%
%     - vol_size:
%       The size of fmri.vol (NIFTI) or fmri.dtseries (CIFTI) or FC.
%
% Example:
% [~, vol, ~] = CBIG_MSHBM_read_fmri('lh.subj01.fsaverage5_roifsaverage3_BI.surf2surf_profile_cen.nii.gz');
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if (~isempty(strfind(fmri_name, 'profile.mat')))
    load(fmri_name);
    vol = single(profile_mat);
    vol_size = size(vol);
    profile_mat = [];
    fmri = [];
elseif (isempty(strfind(fmri_name, '.dtseries.nii')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = single(fmri.vol);
    vol_size = size(vol);
    vol = reshape(vol, prod(vol_size(1:3)), prod(vol_size)/prod(vol_size(1:3)));
    fmri.vol = [];
else
    % if input file is CIFTI file
    fmri = ft_read_cifti(fmri_name);
    vol = single(fmri.dtseries);
    vol_size = size(vol);
    fmri.dtseries = [];
end


end