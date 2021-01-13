function [fmri, vol, vol_size] = CBIG_SPGrad_read_fmri(fmri_name)

% [fmri, vol, vol_size] = CBIG_SPGrad_read_fmri(fmri_name)
% Given the name of functional MRI file (fmri_name), this function read in
% the fmri structure and the content of signals (vol).
% 
% Input:
%     - fmri_name:
%       The full path of input file name.
%
% Output:
%     - fmri:
%       The structure read in by MRIread() or ft_read_cifti(). To save
%       the memory, fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) is
%       set to be empty after it is transfered to "vol".
%
%     - vol:
%       A num_voxels x num_timepoints matrix which is the content of
%       fmri.vol (for NIFTI) or fmri.dtseries (for CIFTI) after reshape.
%
%     - vol_size:
%       The size of fmri.vol (NIFTI) or fmri.dtseries (CIFTI).
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if (isempty(strfind(fmri_name, '.dtseries.nii')))
    % if input file is NIFTI file
    fmri = MRIread(fmri_name);
    vol = fmri.vol;
    vol_size = size(vol);
    if(length(vol_size) == 3)
        vol = reshape(vol, prod(vol_size(1:3)), 1);
    end
    if(length(vol_size) == 4)
        vol = reshape(vol, prod(vol_size(1:3)), vol_size(4));
    end
    fmri.vol = [];
else
    % if input file is CIFTI file
    fmri = ft_read_cifti(fmri_name);
    vol = fmri.dtseries;
    vol = vol(fmri.brainstructure==1|fmri.brainstructure==2, :);
    vol_size = size(vol);
    fmri.dtseries = [];
end