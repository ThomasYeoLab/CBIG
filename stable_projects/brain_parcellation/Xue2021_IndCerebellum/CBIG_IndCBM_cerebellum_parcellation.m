function cbm_labels = CBIG_IndCBM_cerebellum_parcellation(surf_label_file, vol2surf_fc, template_file, output_dir, ...
    varargin)

% cbm_labels = CBIG_IndCBM_cerebellum_parcellation(surf_label_file, vol2surf_fc, template_file, output_dir, varargin)
%
% This function create individual cerebellum parcellation based on it's
% individual cerebral cortical parcellation using winner-take-all algorithm.
% For each cerebellar voxel, this function first identifies the X highest 
% correlated surface vertices. Then, the most frequently assigned network 
% among the X surface vertices is chosen as the "winner" network (each of 
% these X surface vertices is already assigned to a network according to your 
% individual-specific cerebral cortical parcellation). Lastly, we assign the 
% "winner" network to the cerebellar voxel.
% The confidence of the assignment for each cerebellar voxel is measured by
% 1 - N_second_frequent_network/N_most_frequent_network
%
% A cifti template file is required. This file can be created by 
% CBIG_IndCBM_create_template.sh.
%
% This function can work on fsaverage5, fsaverage6, fsaverage, fs_LR_32k,
% fs_LR_164k or your customized resolution. The resolution of labels should
% match the surface resolution of the template file. For customizing 
% surface mesh, please see CBIG_IndCBM_create_template.sh for details.
%
% Input:
%
%     - surf_label_file:
%           mat file of the individual-specific cerebral cortical
%           parcellation. This file should include two variables 
%           lh_labels (Nx1) and rh_labels (Nx1).
%           Example:
%           './examples/example_files/Sub1_example_surf_labels.mat'
%
%     - vol2surf_fc:
%           M x N functional connectivity matrix. 
%           M: cerebellar voxels, same order with the cifti template.
%           N: cerebral cortical vertices, lh and rh concatenated.
%           If this matrix is too large and reading the whole matrix takes
%           too much time and memory, we suggest the user save 'vol2surf_fc' 
%           in a .mat file with 'v7.3'. This function will use matlab 
%           function matfile.m to read functional connectivity. Please
%           refer to CBIG_IndCBM_compute_vol2surf_fc.m for details.
%
%     - template_file: 
%           Cifti template dscalar file for a given mesh and specifying the
%           cerebellar voxels. Please create this file before you run this 
%           function. See CBIG_IndCBM_create_template.sh for how to create 
%           this template file. 
%           Example:
%           './examples/example_files/Sub1_fsaverage5_cerebellum_template.dscalar.nii'
%
%     - output_dir: 
%           Output directory to save parcellation results:
%           dlabel file: Parcellation
%           dscalar file: Confidence of the parcellation. Measured by 
%           1 - N_second_frequent_network/N_most_frequent_network
%
% Optional input:
%
%     - topX: (double or string)
%           Number of cerebral cortical vertices considered. We take X most 
%           correlated surface vertices and assign the most frequent network
%           to the voxel. For example, 'topX', 100 or 'topX', '100'. 
%           If not passed in, default is '100'.
%
%     - output_name: 
%           Specify the output file name. 
%           Example: 'output_name', 'Sub1_parcellation'. 
%
%     - mask_file: 
%           If you need the parcellation in a nifti format, please pass in
%           the nifti volume mask file used in creating cifti template. 
%           Then the function will also generate a nifti file of the 
%           cerebellar parcellation.
%
%     - confidence_off: (double or string)
%           'confidence_off', '1' or 'confidence_off', 1 means confidence 
%           will not be saved. Default is 0.
%
% Output:
%
%     - cbm_labels: 
%           M x 1 vector. Parcellation labels of the cerebellum.
%
% Example:
% CBIG_IndCBM_cerebellum_parcellation('proj/surf_labels.mat', 'proj/sub1_fc', 'proj/template.dscalar.nii', ...
% 'proj/parcellation', 'topX', '200', 'output_name', 'Sub1_parcellation', 'confidence_off', '1')
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Set default optional parameters
pnames = {'topX' 'output_name' 'mask_file' 'confidence_off'};
dflts =  {'100' [] [] 0};
[topX, output_name, mask_file, confidence_off] = internal.stats.parseArgs(pnames, dflts, varargin{:});
if(isempty(output_name))
    output_name = fullfile(output_dir, ['IndCBM_parcellation_top' topX]);
else
    output_name = fullfile(output_dir, output_name);
end

% Change topX to char and confidence_off to num
if(~ischar(topX))
    topX = num2str(topX);
end
if(ischar(confidence_off))
    confidence_off = str2num(confidence_off);
end

addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', ...
'Xue2021_IndCerebellum', 'lib')));

% Load cerebral cortical labels
surf = load(surf_label_file);
if(size(surf.lh_labels, 2) > 1)
    error('lh_labels should be a column vector.');
end
if(size(surf.rh_labels, 2) > 1)
    error('rh_labels should be a column vector.');
end
surf_labels = [surf.lh_labels; surf.rh_labels];

% Run winner-take-all
disp(['Computing cerebellar parcellation with top ' topX ' functionally connected vertices...']);
[cbm_labels, confidence] = CBIG_IndCBM_wta(surf_labels, vol2surf_fc, topX);

% Save parcellation results
disp('Parcellation generated successfully. Saving parcellation results...');
if(~exist(output_dir))
    mkdir(output_dir);
end

% Generate color table. Use 'colors' in surf_label_file if included. Else
% generate a random color table. 
k = length(unique(surf_labels)) - 1;
if(~isfield(surf, 'colors'))
    colors = randi([0, 255], k, 3);
    disp('Generate random color table.');
else
    if(size(surf.colors) == k)
        colors = surf.colors;
    elseif(size(surf.colors, 1) == (k+1))
        disp('We assume the first row of the colortable is medial wall.');
        colors = surf.colors(2:end, :);
    else
        disp('Colortable length does not match label number. Generate random color table instead.');
        colors = randi([0, 255], k, 3);
    end
end
colortable.table = colors;
colortable.struct_names = cell(k, 1);
for i = 1:k
    colortable.struct_names{i} = ['Network_' num2str(i)];
end

% Write dlabel.nii file
CBIG_IndCBM_write_cerebellum_dlabel(surf_labels, cbm_labels, template_file, colortable, output_name);
cifti_file = [output_name '.dlabel.nii'];
disp(['Cifti parcellation saved as ' cifti_file]);

% Write nii.gz file
if(~isempty(mask_file))
    nifti_file = [output_name, '.nii.gz'];
    CBIG_IndCBM_cifti2nifti(cifti_file, mask_file, nifti_file, colortable);
    disp(['Nifti parcellation saved as ' nifti_file]);
end

% Write confidence information
if(~confidence_off)
    cifti = ft_read_cifti(template_file);
    cifti.dscalar = double([zeros(size(surf_labels)); confidence]);
    con_name = [output_name '_confidence'];
    ft_write_cifti(con_name, cifti, 'parameter', 'dscalar');
    con_file = [con_name '.dscalar.nii'];
    disp(['Cifti confidence saved as ' con_file]);
    if(~isempty(mask_file))
        con_nifti_file = [output_name '_confidence.nii.gz']; 
        CBIG_IndCBM_cifti2nifti(con_file, mask_file, con_nifti_file);
        disp(['Nifti confidence saved as ' con_nifti_file]);
    end
end

rmpath(genpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', ...
'Xue2021_IndCerebellum', 'lib')));

end
