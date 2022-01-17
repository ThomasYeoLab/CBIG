function CBIG_ArealMSHBM_generate_radius_mask_Schaefer(resolution, mesh, out_dir)

% CBIG_ArealMSHBM_generate_radius_mask_Schaefer(resolution, mesh, out_dir)
%
% This script will generate radius mask for Schaefer2018 parcellations in a given
% surface space <mesh>.
%
% Input:
%   - resolution: (string)
%     The resolution of the Schaefer2018 parcellation. For example, '400'.
%
%   - mesh: (string)
%     The mesh name of the parcellation surface space. For example, 'fsaverage6', 'fs_LR_32k'.
%
%   - out_dir: (string)
%     The spatial radius mask sparse variables <lh_boundary> and <rh_boundary>
%     will be saved in:
%     <out_dir>/spatial_mask/spatial_mask_<mesh>.mat
%
% Example:
%   CBIG_ArealMSHBM_generate_radius_mask_Schaefer('100', 'fs_LR_32k', '/path_of_out_dir');
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(~isempty(strfind(mesh, 'fsaverage')))
    lh_group_labels = CBIG_read_annotation(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
    'brain_parcellation', 'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3',... 
    mesh, 'label', ['lh.Schaefer2018_' resolution 'Parcels_Kong2022_17Networks_order.annot']));
    rh_group_labels = CBIG_read_annotation(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
    'brain_parcellation', 'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3',... 
    mesh, 'label', ['rh.Schaefer2018_' resolution 'Parcels_Kong2022_17Networks_order.annot']));
    lh_labels = lh_group_labels - 1;
    rh_labels = rh_group_labels - 1;
elseif(~isempty(strfind(mesh, 'fs_LR')))
    group_labels = ft_read_cifti(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects',...
    'brain_parcellation', 'Schaefer2018_LocalGlobal', 'Parcellations', 'HCP', 'fslr32k', 'cifti',...
      ['Schaefer2018_' resolution 'Parcels_Kong2022_17Networks_order.dscalar.nii']));
    lh_labels = group_labels.dscalar(1:32492);
    rh_labels = group_labels.dscalar(32493:end);
    rh_labels(rh_labels~=0) = rh_labels(rh_labels~=0) - max(lh_labels);
end
CBIG_ArealMSHBM_generate_radius_mask(lh_labels, rh_labels, mesh, '30', out_dir);

end



