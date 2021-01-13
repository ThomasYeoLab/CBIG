function CBIG_IndCBM_example_wrapper(output_dir)

% CBIG_IndCBM_example_wrapper(output_dir)
% This function use fake data to create individual cerebellum parcellation
% based on a group-level cerebral cortical parcellation.
%
% Input:
%     - output_dir: 
%           Results will be saved under this folder. Final cerebellum
%           parcellation will be under output_dir/parcellation
%
% Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% Add code path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
IndCBM_code_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Xue2021_IndCerebellum');
IndCBM_example_dir = fullfile(IndCBM_code_dir, 'examples');
MSHBM_code_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM');
MSHBM_output_dir = fullfile(output_dir, 'MSHBM');
addpath(IndCBM_code_dir);
addpath(genpath(MSHBM_code_dir));

% Remove the output folder if exist
if(exist(output_dir, 'dir'))
    disp('Output folder already exists. Remove output folder.');
    rmdir(output_dir, 's')
end

% Create data list for example under output_dir/list and output_dir/MSHBM
cmd = ['sh ' fullfile(IndCBM_example_dir, 'CBIG_IndCBM_generate_example_list.sh') ' ' output_dir];
system(cmd);

% Cerebral cortical parcellation
% Compute surface profile on fsaverage5 mesh
mesh = 'fsaverage5';
cmd = ['sh ' fullfile(IndCBM_code_dir, 'CBIG_IndCBM_compute_profile.sh') ' ' mesh ' ' MSHBM_output_dir];
system(cmd);

% Average surface profile for 2 subjects
lh_profile = fullfile(MSHBM_output_dir, 'profiles', 'avg_profile', 'lh_fsaverage5_roifsaverage3_avg_profile.nii.gz');
rh_profile = fullfile(MSHBM_output_dir, 'profiles', 'avg_profile', 'rh_fsaverage5_roifsaverage3_avg_profile.nii.gz');

% Generate MSHBM parameters using average profile
load(fullfile(IndCBM_example_dir, 'group_surf', 'fsaverage5_10networks'));
clustered = CBIG_IndCBM_generate_MSHBM_params(lh_profile, rh_profile, lh_labels, rh_labels);

% Save parameters to output_dir/MSHBM/group
output_group_dir = fullfile(MSHBM_output_dir, 'group');
if(~exist(output_group_dir))
    mkdir(output_group_dir);
end
if(exist('colors', 'var'))
    save(fullfile(MSHBM_output_dir, 'group', 'group.mat'), 'lh_labels', 'rh_labels', 'colors', 'clustered');
else
    save(fullfile(MSHBM_output_dir, 'group', 'group.mat'), 'lh_labels', 'rh_labels', 'clustered');
end

% Generate individual surface parcellation using MSHBM algorithm. Results
% will be saved under output_dir/MSHBM/ind_parcellation
CBIG_MSHBM_estimate_group_priors(MSHBM_output_dir, 'fsaverage5', '2', '2', '10', 'conv_th', '1e-5', ...
    'save_all', '1', 'max_iter', '5'); 
CBIG_IndCBM_extract_MSHBM_result(MSHBM_output_dir);

% Cerebellar parcellation
for i = 1:2
    sub = ['sub' num2str(i)];
    
    % Create cifti template for individual subjects based on a cerebellum 
    % mask in the volume. Template file will be saved under output_dir/template
    mask_file = fullfile(IndCBM_example_dir, 'input', 'mask', sub, [sub '_bin_mask_4mm.nii.gz']);
    cmd = ['sh ' fullfile(IndCBM_code_dir, 'CBIG_IndCBM_create_template.sh') ' ' mesh ' ' mask_file ' '];
    cmd = [cmd fullfile(output_dir, 'template') ' -m 23'];
    system(cmd);

    % Compute functional connectivity between cerebellum and cerebral cortex
    lh_surf_file = fullfile(output_dir, 'list', ['lh_' sub '_list.txt']);
    rh_surf_file = fullfile(output_dir, 'list', ['rh_' sub '_list.txt']);
    vol_file = fullfile(output_dir, 'list', ['vol_' sub '_list.txt']);
    template_file = fullfile(output_dir, 'template', [sub '_bin_mask_4mm_fsaverage5_cerebellum_template.dscalar.nii']);
    vol2surf_fc = CBIG_IndCBM_compute_vol2surf_fc(lh_surf_file, rh_surf_file, vol_file, template_file);
    
    % Use functional connectivity and individual surface parcellation to
    % create cerebellar parcellation
    out_p_dir = fullfile(output_dir, 'parcellation', sub);
    if(~exist(out_p_dir))
        mkdir(out_p_dir);
    end
    surf_label_file = fullfile(MSHBM_output_dir, 'ind_parcellation', ['Ind_parcellation_MSHBM_' sub]);
    CBIG_IndCBM_cerebellum_parcellation(surf_label_file, vol2surf_fc, template_file, out_p_dir, 'mask_file', mask_file);
end

disp('Example results generated successfully.');

% Remove path
rmpath(IndCBM_code_dir);
rmpath(genpath(MSHBM_code_dir));

end
