function CBIG_SPGrad_construct_cifti_format(mesh, mask, data, filename)

% CBIG_SPGrad_construct_cifti_format(mesh, mask, data, filename)
%
% This script write data as cifti format
%
% Input:
%     - mesh:
%       resolution of surface mesh, e.g. 'fsaverage6', 'fs_LR_32k'
%
%     - medial_mask: (#num_vertices x 1 binary vector or 'NONE')
%       the medial area mask. A #num_vertices x 1 binary vector. The medial area should be denoted 
%       as 1, others should be denoted as 0.
%
%     - data: (#num_vertices x K)
%       data matrix that will be saved as cifti format
%
%     - filename:
%       output filename. For example, '/path/your_data.dtseries.nii'.
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if(strcmp(mesh,'fs_LR_32k'))
    cifti_template = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'speedup_gradients', ...
    'utilities', 'fslr_surface_template', 'cifti_template.dscalar.nii');
    cifti = ciftiopen(cifti_template, 'wb_command');
    cifti.cdata = [];
    cifti.cdata(~mask,:) = data;
    ciftisavereset(cifti, filename, 'wb_command');
    tmp = ft_read_cifti_mod(filename);
    tmp.data(mask,:) = [];
    tmp.brainstructure(mask) = -1;
    ft_write_cifti_mod(filename,tmp);

elseif(~isempty(strfind(mesh,'fsaverage6')))
    cifti_template = fullfile(getenv('CBIG_CODE_DIR'), 'utilities', 'matlab', 'speedup_gradients', ...
    'utilities', 'fs6_surface_template', 'cifti_template.dscalar.nii');
    cifti = ciftiopen(cifti_template, 'wb_command');
    cifti.cdata = [];
    cifti.cdata(~mask,:) = data;
    ciftisavereset(cifti, filename, 'wb_command');
    tmp = ft_read_cifti_mod(filename);
    tmp.data(mask,:) = [];
    tmp.brainstructure(mask) = -1;
    ft_write_cifti_mod(filename,tmp);
else
    error(['We do not have midthickness surface mesh for ' mesh 'in current version.']);
end

