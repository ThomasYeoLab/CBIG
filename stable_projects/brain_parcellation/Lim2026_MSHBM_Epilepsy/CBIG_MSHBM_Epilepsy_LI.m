function [lh_labels, rh_labels, LI_lang] = CBIG_MSHBM_Epilepsy_LI(params)
% ====================================
%% Help
% ====================================
%
% CBIG_MSHBM_Epilepsy_LI
%
% Estimates an individual-specific 15-network cortical parcellation for a
% patient with drug-resistant epilepsy using the NIH 34-subject MS-HBM
% model, and computes the language laterality index (LI) from the
% individual parcellation. Requires Kong2019_MSHBM scripts as dependencies.
%
% The individual parcellation is saved to:
%   params.project_dir/ind_parcellation/test_set/
%     Ind_parcellation_MSHBM_sub1_w80_MRF10.mat
% The LI is saved to:
%   params.project_dir/LI.mat (variable: LI_lang)
%
% INPUTS:
%   params: A structure with the following fields:
%     - project_dir (required):
%       String. The output directory for this subject. Must be unique per
%       subject. The individual parcellation and LI.mat are saved here.
%
%     - lh_fMRI_list (required):
%       1xN or Nx1 cell array of strings. Full paths to the left-hemisphere
%       fMRI surface files for each session. Can be a string for one
%       session.
%
%     - rh_fMRI_list (required for fsaverage6/fsaverage5):
%       1xN or Nx1 cell array of strings. Full paths to the right-hemisphere
%       fMRI surface files for each session. Can be a string for one
%       session.
%
%     - censor_list (required):
%       1xN or Nx1 cell array of strings. Full paths to the motion censor
%       files for each session (one row per time point: 1=keep, 0=censor).
%       Set to 'NONE' if no censoring is required. Can be a string for one
%       session.
%
%     - target_mesh (optional):
%       String. Surface mesh. Default: 'fsaverage6'.
%
%     - group_prior (optional):
%       String. Full path to the group prior Params_Final.mat. Defaults to
%       the NIH 34-subject drug-resistant epilepsy MS-HBM model:
%         Lim2026_MSHBM_Epilepsy/MSHBM_Epilepsy/Params_Final.mat
%
%     - w (optional):
%       String. Weight for the group spatial prior. Default: '80'.
%
%     - c (optional):
%       String. Weight for the MRF smoothness constraint. Default: '10'.
%
%     - overwrite_flag (optional):
%       0 or 1. If 1, deletes project_dir before running. Default: 0.
%
% OUTPUTS:
%   lh_labels: Vector. Network labels for left hemisphere vertices.
%   rh_labels: Vector. Network labels for right hemisphere vertices.
%   LI_lang:   Scalar. Language laterality index:
%                (LH language vertices - RH language vertices) /
%                (LH language vertices + RH language vertices)
%              Language networks in the NIH model:
%                Lang A (primary)   = label 8
%                Lang B (secondary) = label 11
%              Positive LI = left-lateralised; negative = right-lateralised.
%
% EXAMPLE:
%   params.project_dir  = '/myproject/patient01';
%   params.lh_fMRI_list = {'/mydata/lh.sess1.nii.gz', '/mydata/lh.sess2.nii.gz'};
%   params.rh_fMRI_list = {'/mydata/rh.sess1.nii.gz', '/mydata/rh.sess2.nii.gz'};
%   params.censor_list  = {'/mydata/censor1.txt', '/mydata/censor2.txt'};
%   [lh_labels, rh_labels, LI_lang] = CBIG_MSHBM_Epilepsy_LI(params);
%
% Written by Mervyn Lim Jun Rui and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% ====================================
%% Setup
% ====================================

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
lim2026_dir = fileparts(mfilename('fullpath'));

addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM'));

% Default group prior: NIH 34-subject epilepsy MS-HBM model
if ~isfield(params, 'group_prior')
    params.group_prior = fullfile(lim2026_dir, 'MSHBM_Epilepsy', 'Params_Final.mat');
end

% Epilepsy-optimised parcellation parameters
if ~isfield(params, 'target_mesh')
    params.target_mesh = 'fsaverage6';
end
if ~isfield(params, 'w')
    params.w = '80';
end
if ~isfield(params, 'c')
    params.c = '10';
end

% ====================================
%% Estimate Individual Parcellation
% ====================================

[lh_labels, rh_labels] = CBIG_MSHBM_parcellation_single_subject(params);

% ====================================
%% Compute Language Laterality Index
% ====================================

% Language network labels for the NIH 34-subject MS-HBM epilepsy model:
%   Lang A (primary language network)   = label 8
%   Lang B (secondary language network) = label 11
lh_lang = sum(lh_labels == 8 | lh_labels == 11);
rh_lang = sum(rh_labels == 8 | rh_labels == 11);
LI_lang = (lh_lang - rh_lang) / (lh_lang + rh_lang);

fprintf('Language LI (Lang A + B): %.4f\n', LI_lang);
save(fullfile(params.project_dir, 'LI.mat'), 'LI_lang');
fprintf('LI saved to %s\n', fullfile(params.project_dir, 'LI.mat'));

% ====================================
%% Cleanup
% ====================================

rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM'));

end
