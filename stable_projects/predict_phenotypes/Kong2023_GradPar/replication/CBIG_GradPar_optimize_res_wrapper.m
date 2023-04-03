function CBIG_GradPar_optimize_res_wrapper(LRR_out_dir, KRR_out_dir, opt_res_out_dir)

% CBIG_GradPar_optimize_res_wrapper(LRR_out_dir, KRR_out_dir)
%
% This script is a wrapper function to optimize the resolution for sICA and Kong2021. Before running this
% script, please make sure you have run the following scripts to generate prediction results for sICA and
% Kong2021 across different resolutions:
% sh ./CBIG_GradPar_LRR_frac_each_resolution_submit.sh <LRR_out_dir>
% sh ./CBIG_GradPar_KRR_each_resolution_submit.sh <KRR_out_dir>
% Please note that this wrapper is designed to replicate a subset of results of Kong2023 to reduce time.
% The script can be modified to replicate all results.
%
% Input:
%   - LRR_out_dir:
%       The output directory of LRR fracridge prediction results by running
%       sh ./CBIG_GradPar_LRR_frac_each_resolution_submit.sh <LRR_out_dir>
%
%   - KRR_out_dir:
%       The output directory of KRR prediction results by running
%       sh ./CBIG_GradPar_KRR_each_resolution_submit.sh <KRR_out_dir>
%
%   - opt_res_out_dir:
%       The output directory to save the optimized resolution results.
%
% Output:
%   - The prediction resilts with optimal resolutions for sICA using LRR is under: 
%     <opt_res_out_dir>/sICA_opt_res
%   - The prediction resilts with optimal resolutions for Kong2021 using KRR is under:
%     <opt_res_out_dir>/Kong2021_opt_res
%
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar'));
REPDATA_DIR = fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar');

% Set up project_set
% Please note here to reduce the running time, we only use:
% 2 resolutions for sICA in ABCD dataset using LRR
% 2 resolutions for Kong2021 in HCP dataset using KRR.
% The projects paths and input parameters were set up based on the set up in:
% CBIG_GradPar_LRR_frac_each_resolution_submit.sh
% CBIG_GradPar_KRR_each_resolution_submit.sh
% To replicate all results, please generate prediction results for all parcellations/gradient approaches with
% all resolutions for both KRR and LRR in HCP and ABCD dataset. And change the projects set to include paths.
LRR_projects_set_ABCD = {fullfile(LRR_out_dir, 'ABCD', 'sICA', '50'), ...
    fullfile(LRR_out_dir, 'ABCD', 'sICA', '100')};
KRR_projects_set_HCP = {fullfile(KRR_out_dir, 'HCP', 'Kong2021', '100'), ...
    fullfile(KRR_out_dir, 'HCP', 'Kong2021', '200')};

% Set up parameters for ABCD, sICA, LRR
params_ABCD.num_splits = '';
params_ABCD.num_folds = '84';
params_ABCD.num_behaviors = '39';
params_ABCD.outstem = '36_behaviors_3_components';
params_ABCD.project_set = LRR_projects_set_ABCD;
params_ABCD.out_dir = fullfile(opt_res_out_dir, 'sICA_opt_res');

% Optimize resolution
CBIG_GradPar_LRR_frac_optimize_res_wrapper(params_ABCD);

% Set up parameters for HCP, Kong2021, KRR
params_HCP.nsplit = '1';
params_HCP.num_folds = '20';
params_HCP.num_behaviors = '61';
params_HCP.outstem = '58_behaviors_3_components';
params_HCP.project_set = KRR_projects_set_HCP;
params_HCP.out_dir = fullfile(opt_res_out_dir, 'Kong2021_opt_res');
lambda_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params', 'lambda_set.mat'));
params_HCP.lambda_set = lambda_path.lambda_set;
ker_param_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params', 'ker_param.mat'));
params_HCP.ker_param = ker_param_path.ker_param;
sub_fold_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params',...
    params_HCP.nsplit, 'no_relative_S1200_20_fold_sub_list.mat'));
params_HCP.sub_fold = sub_fold_path.sub_fold;

CBIG_GradPar_KRR_optimize_res_wrapper(params_HCP);

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar'));

end