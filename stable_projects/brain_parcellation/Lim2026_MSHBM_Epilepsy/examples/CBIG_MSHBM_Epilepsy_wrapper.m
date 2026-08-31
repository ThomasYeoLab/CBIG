function CBIG_MSHBM_Epilepsy_wrapper(out_dir)

% ====================================
%% Help
% ====================================
%
% CBIG_MSHBM_Epilepsy_wrapper
%
% Generates an individual-specific 15-network cortical parcellation for a
% single epilepsy patient using the NIH 34-subject drug-resistant epilepsy
% MS-HBM model, then computes the language laterality index (LI) from the
% parcellation. This wrapper demonstrates how to call
% CBIG_MSHBM_Epilepsy_LI using the 5 example resting-state fMRI runs and
% corresponding motion censor files provided with the code release.
%
% INPUTS:
%   out_dir (optional): Full path to the output directory for this subject.
%   Must be unique per subject. Intermediate files (FC profiles, data
%   lists), the final parcellation, and the LI are saved here. Defaults to
%   out_dir/ in the same folder as this wrapper script (i.e.
%   Lim2026_MSHBM_Epilepsy/examples/out_dir/).
%   To check outputs against the reference, use
%   CBIG_MSHBM_Epilepsy_check_example_results(out_dir).
%
% OUTPUTS:
%   Individual parcellation labels:
%     out_dir/ind_parcellation/test_set/
%       Ind_parcellation_MSHBM_sub1_w80_MRF10.mat
%   Brain surface visualisation:
%     out_dir/epilepsy_w80c10.png
%   Language laterality index:
%     out_dir/LI.mat  (variable: LI_lang)
%
% EXAMPLE USAGE:
%   CBIG_MSHBM_Epilepsy_wrapper();
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Prepare Inputs
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');

% Locate directories relative to this wrapper's own location
examples_dir = fileparts(mfilename('fullpath'));
lim2026_dir  = fileparts(examples_dir);

% Default output directory
if nargin < 1 || isempty(out_dir)
    out_dir = fullfile(examples_dir, 'out_dir');
end

% Output directory for this subject
params.project_dir = out_dir;

% Group prior: NIH 34-subject drug-resistant epilepsy MS-HBM model
params.group_prior = fullfile(lim2026_dir, 'MSHBM_Epilepsy', 'Params_Final.mat');

colors_file = fullfile(lim2026_dir, 'MSHBM_Epilepsy', 'NIH_34_group', 'group_colors.mat');

% fMRI surface files: 5 example rs-fMRI runs
lh_files = dir(fullfile(examples_dir, 'input', 'lh.example_bld*_fs6_sm6.nii.gz'));
params.lh_fMRI_list = {};
for i = 1:length(lh_files)
    params.lh_fMRI_list{i} = fullfile(lh_files(i).folder, lh_files(i).name);
end

rh_files = dir(fullfile(examples_dir, 'input', 'rh.example_bld*_fs6_sm6.nii.gz'));
params.rh_fMRI_list = {};
for i = 1:length(rh_files)
    params.rh_fMRI_list{i} = fullfile(rh_files(i).folder, rh_files(i).name);
end

% Censor files: binary vectors (1=keep, 0=censor) for each run
censor_files = dir(fullfile(examples_dir, 'input', 'example_bld*_motion_outliers.txt'));
params.censor_list = {};
for i = 1:length(censor_files)
    params.censor_list{i} = fullfile(censor_files(i).folder, censor_files(i).name);
end

% Surface mesh and parcellation parameters
params.target_mesh = 'fsaverage6';
params.w = '80';
params.c = '10';

% ====================================
%% Run Parcellation and Compute LI
% ====================================

addpath(lim2026_dir);
[lh_labels, rh_labels, LI_lang] = CBIG_MSHBM_Epilepsy_LI(params);
rmpath(lim2026_dir);

% ====================================
%% Save Figure
% ====================================

addpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab', 'figure_utilities'));
group = load(colors_file);
CBIG_DrawSurfaceMaps(lh_labels, rh_labels, 'fsaverage6', 'inflated', 0, 15, group.colors);
saveas(gcf, fullfile(out_dir, ['epilepsy_w' params.w 'c' params.c '.png']));
close(gcf);
rmpath(fullfile(CBIG_CODE_DIR, 'utilities', 'matlab', 'figure_utilities'));

fprintf('Successfully saved parcellation figure for example subject.\n');
