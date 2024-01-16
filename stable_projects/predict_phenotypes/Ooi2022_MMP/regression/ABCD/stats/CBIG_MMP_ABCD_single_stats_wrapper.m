function p_value = CBIG_MMP_ABCD_single_stats_wrapper(outdir, stats_dir)

% p_value = CBIG_ABCD_single_stats_wrapper(outdir, stats_dir)
%
% Wrapper function that calculates whether each model performs better than chance.
% Assumes that all regression outputs are in the same directory.
%
% Inputs:
%
%   - outdir
%     Directory to results of regression models.
%
%   - outdir
%     Directory to permutations of null models
%
% Outputs:
%   - p-value
%     A matrix of #behav_ind x #models. P-value for each model indicating in predicting
%     specified behaviour.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% set utility directory
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'utilities'))

% list of models
models = {'KRR_features_cv' 'KRR_features_ca' 'KRR_features_ct' 'KRR_features_tbss_FA' ...
    'KRR_features_tbss_MD' 'KRR_features_tbss_AD' 'KRR_features_tbss_RD' 'KRR_features_tbss_OD' ...
    'KRR_features_tbss_ICVF' 'KRR_features_tbss_ISOVF' 'KRR_features_schaefer_FA' 'KRR_features_schaefer_MD' ...
    'KRR_features_schaefer_AD' 'KRR_features_schaefer_RD' 'KRR_features_schaefer_OD' ...
    'KRR_features_schaefer_ICVF' 'KRR_features_schaefer_ISOVF' 'KRR_features_schaefer_streamcount_log' ...
    'KRR_features_schaefer_streamlen'  'KRR_features_rs' 'KRR_features_mid' ...
    'KRR_features_sst' 'KRR_features_nback' 'multiKRR_k1_rs_k2_mid_k3_sst_k4_nback' ...
    'stacking_LRR_rs_mid_sst_nback' 'stacking_LRR_all'};
metric = 'corr';
N_folds = 120;
behav_ind = [37:39];
N_behav = length(behav_ind);
num_perms = 10000;

% save p_values
clear p_value
for i = 1:length(models)
    outstem = models{i};
    ref_dir = fullfile(outdir,outstem);
    perm_dir = fullfile(stats_dir, outstem);
    p_value(i,:) = CBIG_MMP_compute_ABCD_permutation_p_value(outdir, behav_ind, outstem, perm_dir, num_perms);
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'utilities'))
