function [p_single, pairwise_stats, threshold] = CBIG_MMP_ABCD_FDR_wrapper

% ABCD wrapper for calculating p-values of whether models perform above chance, and the pairwise comparison of
% models. FDR is performed for all p-values subsequently. Outputs final threshold for statistics.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% define output and stats directories: MODIFY FOR YOUR OWN OUTPUT DIRECTORY
out_dir = fullfile('/home', 'leon_ooi', 'storage', 'Multimodal_prediction_project', ...
    'replication', 'ABCD', 'output');
stats_dir = fullfile('/home', 'leon_ooi', 'storage', 'Multimodal_prediction_project', ...
    'replication', 'ABCD', 'stats');

% add stats path
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'stats'))
% individual model stats
p_single = CBIG_MMP_ABCD_single_stats_wrapper(out_dir, stats_dir);
% model comparison stats
[pairwise_stats, p_compbehav_tmp, p_combined_tmp] = CBIG_MMP_ABCD_pairwise_stats_wrapper(out_dir);
p_combined_noCOD = p_combined_tmp(1:5,:);
p_single_noBest = p_compbehav_tmp(1:9,:);
p_rearr = [p_single_noBest(:); p_combined_noCOD(:); p_single(:)];

% perform FDR
[ind, threshold] = FDR(p_rearr, 0.05);

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'ABCD', 'stats'))
