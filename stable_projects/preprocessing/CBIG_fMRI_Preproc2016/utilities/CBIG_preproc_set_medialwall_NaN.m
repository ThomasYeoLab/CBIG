function CBIG_preproc_set_medialwall_NaN( hemi, surf_mesh, input_name, output_name )

% CBIG_preproc_set_medialwall_NaN( hemi, surf_mesh, input_name, output_name )
% 
% To force users do not use the timeseries within medial wall, we will set
% all medial wall timeseries of surface fMRI data (input_name) to be NaN in
% this function, and save out to "output_name". The medial wall is defined
% by $FREESURFERHOME/subjects/<surf_mesh>/label/?h.Medial_wall.label.
% 
% Inputs:
%     - hemi:
%       Hemisphere of input, 'lh' or 'rh'.
% 
%     - surf_mesh:
%       The space resolution, e.g., 'fsaverage6'.
% 
%     - input_name:
%       The full filename of input surface data, 
%       e.g., <path to data>/lh.Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6.nii.gz 
% 
%     - output_name:
%       The full filename of output_surface data, 
%       e.g., <path to data>/lh.Sub0001_Ses1_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_medialwallNaN.nii.gz
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


hemi_medialwall_info = read_label(surf_mesh ,[hemi '.Medial_wall']); 
% hemi_medialwall_info is a 5-column matrix. The first column is the medial
% wall vertices' indices. The index starts from 0. To transform to matlab index,
% it needs to add 1.
medialwall_ind = hemi_medialwall_info(:, 1) + 1;

mri = MRIread(input_name);
size_vol = size(mri.vol);
vol = reshape(mri.vol, size_vol(1)*size_vol(2) * size_vol(3), size_vol(4));

vol(medialwall_ind, :) = nan;

mri.vol = reshape(vol, size_vol);
MRIwrite(mri, output_name);

end