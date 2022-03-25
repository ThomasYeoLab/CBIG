function tr = CBIG_preproc_infer_TR(input_fMRI)

% tr = CBIG_preproc_infer_TR(input_vol)
% 
% This function computes the TR (repetition time) of fMRI data
% 
% Inputs:
%     - input_vol:
%       Full path of the input fMRI data
%
% Outputs:
%     - tr
%       TR of the input fMRI in seconds
% 
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
if ~contains(input_fMRI, '.dtseries.nii')
    % if input_vol is a nifti file
    mri = MRIread(input_fMRI,1);
    tr = mri.tr/1000;
else
    % if input_vol is a cifti file
    mri = ft_read_cifti(input_fMRI,'readdata',false);
    tr = mri.time(2) - mri.time(1);
end

end