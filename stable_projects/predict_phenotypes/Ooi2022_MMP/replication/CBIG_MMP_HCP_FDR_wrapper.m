function [p_single, pairwise_stats, threshold] = CBIG_MMP_HCP_FDR_wrapper

% HCP wrapper for calculating p-values of whether models perform above chance, and the pairwise comparison of
% models. FDR is performed for all p-values subsequently. Outputs final threshold for statistics.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% define output and stats directories: MODIFY FOR YOUR OWN OUTPUT DIRECTORY
out_dir = fullfile('/home', 'leon_ooi', 'storage', 'Multimodal_prediction_project', ...
    'replication', 'HCP', 'output');
stats_dir = fullfile('/home', 'leon_ooi', 'storage', 'Multimodal_prediction_project', ...
    'replication', 'HCP', 'stats');

% add stats path
addpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'HCP', 'stats'))
% individual model stats
p_single = CBIG_MMP_HCP_single_stats_wrapper(out_dir, stats_dir);
% model comparison stats
[pairwise_stats, p_compbehav_tmp, p_combined_tmp] = CBIG_MMP_HCP_pairwise_stats_wrapper(out_dir);
p_combined_noCOD = p_combined_tmp(1:5,:);
p_single_noBest = p_compbehav_tmp(1:9,:);
p_rearr = [p_single_noBest(:); p_combined_noCOD(:); p_single(:)];

% perform FDR
[ind, threshold] = FDR(p_rearr, 0.05);

rmpath(fullfile(getenv('CBIG_CODE_DIR'),'stable_projects', 'predict_phenotypes', ...
    'Ooi2022_MMP', 'regression', 'HCP', 'stats'))
