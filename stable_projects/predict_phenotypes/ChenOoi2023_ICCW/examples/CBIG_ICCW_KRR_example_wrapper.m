function CBIG_ICCW_KRR_example_wrapper(out_dir)

% CBIG_ICCW_KRR_example_wrapper(out_dir)
%
% This function generates the KRR results for the example in ChenOoi2023.
% Subsequently the Haufe-transformed regression weights will be generated
% and the ICC pf the Haufe-transformed regression weights is calculated.
%
% Input:
%  - out_dir
%    A path where the results of the example will be saved.
%
% Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% define directories and set up dependencies
project_code_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
    'predict_phenotypes', 'ChenOoi2023_ICCW', 'analysis', 'utilities');
parent_dir = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', ...
    'ChenOoi2023_ICCW', 'examples');
input_dir = fullfile(parent_dir, 'input');
addpath(genpath(project_code_dir));

%% STEP 1: run prediction portion
% load input
load(fullfile(input_dir, 'no_relative_2_fold_sub_list.mat'))
load(fullfile(input_dir, 'y.mat'))
load(fullfile(input_dir, 'covariates.mat'))
load(fullfile(input_dir, 'RSFC.mat'))

% prepare param for KRR
param.sub_fold = sub_fold;
param.y = y;
param.feature_mat = corr_mat;
param.covariates = covariates;
param.num_inner_folds = 5;
param.outdir = out_dir;
param.outstem = '2cog';
param.with_bias = 1;
param.ker_param.type = 'corr';
param.ker_param.scale = nan;
param.lambda_set = [0 0.01 0.1 1 10 100 1000];
param.threshold_set = nan;
param.metric = 'corr';
param.cov_X = [];

% create outdir
if(exist(param.outdir, 'dir'))
    rmdir(param.outdir, 's')
end
mkdir(param.outdir)
save(fullfile(param.outdir, 'setup.mat'), '-struct', 'param')

% call the KRR workflow function
CBIG_KRR_workflow_LITE( fullfile(param.outdir, 'setup.mat') )

%% STEP 2: run interpretation
% load results
load(fullfile(param.outdir, 'final_result_2cog.mat'))

% extract lower triangle of FC
FC_vec = zeros(780,50);
lt = logical(tril(ones(40,40),-1));
for subs = 1:size(corr_mat,3)
    FC_tmp = corr_mat(:,:,subs);
    FC_vec(:,subs) = FC_tmp(lt);
end
% calculate Haufe transformed regression weights
pfm1 = CBIG_ICCW_cov_matrix(FC_vec(:,sub_fold(1).fold_index==0)',y_pred_train{1});
pfm2 = CBIG_ICCW_cov_matrix(FC_vec(:,sub_fold(2).fold_index==0)',y_pred_train{2});

% get icc
icc = CBIG_ICCW_ICC_1_1([pfm1 pfm2]);
save(fullfile(param.outdir, 'pfm_icc.mat'), 'pfm1', 'pfm2', 'icc');

rmpath(genpath(project_code_dir));
end
