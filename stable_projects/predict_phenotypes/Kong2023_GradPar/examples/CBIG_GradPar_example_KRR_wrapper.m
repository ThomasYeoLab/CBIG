function CBIG_GradPar_example_LRR_frac_wrapper(out_dir)

% CBIG_GradPar_example_LRR_frac_wrapper(out_dir)
%
% This scripts is a wrapper function to run the LRR_fracfridge workflow using example
% data for two resolutions 100 and 200. 
%
% Inputs:
%   - out_dir
%     The output directory where the results will be saved.
%
%  output (saved)
%   - acc_metric_train - training accuracy
%   - acc_corr_test - test accuracies (given in correlation)
%   - y_predict - predicted target values
%   - optimal_statistics - a cell array of size equal to the number of
%   folds. Each element in the cell array is a structure storing the
%   accuracies of each possible accuracy metric (eg. corr, MAE, etc).
% 
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
 
addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar'));

% Set up path to example data
DATA_DIR = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes',...
    'Kong2023_GradPar', 'examples', 'input');
out_dir = fullfile(out_dir, 'KRR');

% Loop through two resolutions
for res = [100, 200]
    res = num2str(res);
    % Check if the output file already exists
    if(~exist(fullfile(out_dir, res, 'final_result_example_targets.mat')))
        % Set up parameters
        % Set up cross-validation folds
        sub_fold_path = matfile(fullfile(DATA_DIR, 'no_relative_3_fold_sub_list.mat'));
        params.sub_fold = sub_fold_path.sub_fold;

        % Set up feature matrix
        feature_mat_path = matfile(fullfile(DATA_DIR, ['feature_mat_' res '.mat']));
        feature_names = fieldnames(feature_mat_path);
        params.feature_mat = feature_mat_path.(feature_names{2});

        % Set up covariates matrix
        covariates_path = matfile(fullfile(DATA_DIR, 'regressors.mat'));
        params.covariates = covariates_path.regressors;

        % Set up target matrix
        y_path = matfile(fullfile(DATA_DIR, 'targets.mat'));
        params.y = y_path.y_in;

        % Set up output directory
        params.outdir = fullfile(out_dir, res);

        % Set up output file name
        params.outstem = 'example_targets';

        % Set up parameters for LRR_fracridge: number of innerloop cross-validation folds
        params.num_inner_folds = 3;

        % Set up metric used for tuning hyperparameter
        params.metric = 'predictive_COD';

        % Set up bias flag
        params.with_bias = 1;

        % Set up kernel type
        params.ker_param.type = 'corr';
        params.ker_param.scale = NaN;

        % Set up regularization parameter searching range
        params.lambda_set = [ 0 0.00001 0.0001 0.001 0.004 0.007 0.01 0.04 0.07 0.1 0.4 0.7...
            1 1.5 2 2.5 3 3.5 4 5 10 15 20];

        % Call KRR workflow
        CBIG_KRR_workflow_LITE(params);
    else
        disp('File generated!')
    end
end

% Optimize resolution
KRR_projects_set = {fullfile(out_dir, '100'), ...
                    fullfile(out_dir, '200')};

% Set up parameters for optimizing resolution
params_opt.nsplit = '';
params_opt.num_folds = '3';
params_opt.num_behaviors = '3';
params_opt.outstem = 'example_targets';
params_opt.project_set = KRR_projects_set;
params_opt.out_dir = fullfile(out_dir, 'opt_res');
params_opt.lambda_set = params.lambda_set;
params_opt.ker_param = params.ker_param;
params_opt.sub_fold = params.sub_fold;

CBIG_GradPar_KRR_optimize_res_wrapper(params_opt);

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar'));
end