function CBIG_GradPar_KRR_each_resolution_HCP_wrapper(out_dir, project_name, res, nsplit)

% CBIG_GradPar_KRR_each_resolution_HCP_wrapper(out_dir, project_name, res, nsplit)
%
% This script is a wrapper function to run the KRR workflow in HCP
% dataset for each resolution of different approaches Schaefer2018, Kong2021, sICA,
% LocalGrad, and PrincipalGrad.
%
% Inputs:
%   - out_dir
%     The output directory where the results will be saved.
%
%   - project_name 
%     The name of the project. It can be one of the following:
%     'Schaefer2018', 'Kong2021', 'sICA', 'LocalGrad', and 'PrincipalGrad'.
%
%   - res
%     The resolution of the project. It can be one of the following:
%     -- Schaefer2018, Kong2021:
%     '100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'.       
%     -- sICA:
%     '50','100','200','300'.
%     -- LocalGrad:
%     There is only one resolution for LocalGrad, which is '1'.
%     -- PrincipalGrad:
%     '1','5','10','20','40','60','80','100'.
%
%   - nsplit
%     The number of split for the project. We did 100 splits for each project.
%     <nsplit> can be one of the following:'1','2',...,'100'. In each split,
%     we randomly split the subjects into 20 folds. <nsplit> is the random seed
%     number we used to split the subjects for cross-validation. For example, if
%     <nsplit> is '1', then we will use the first split of the subjects.
%
%  output (saved)
%   - acc_corr_train - training accuracy (given in correlation)
%   - acc_metric_test - test accuracies (given in correlation)
%   - y_predict - predicted target values
%   - optimal_statistics - a cell array of size equal to the number of
%   folds. Each element in the cell array is a structure storing the
%   accuracies of each possible accuracy metric (eg. corr, MAE, etc).
% 
% Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
    
% Set up path to replication data
REPDATA_DIR = fullfile(getenv('CBIG_REPDATA_DIR'), 'stable_projects', 'predict_phenotypes', 'Kong2023_GradPar');

% Check if the output file already exists
if(~exist(fullfile(out_dir, project_name, res, nsplit, 'final_result_58_behaviors_3_components.mat')))
    % Set up parameters
    % Set up cross-validation folds for each split
    sub_fold_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params',... 
        nsplit, 'no_relative_S1200_20_fold_sub_list.mat'));
    params.sub_fold = sub_fold_path.sub_fold;

    % Set up feature matrix
    if(~isempty(strfind(project_name , 'PrincipalGrad')))
        % To save space, we only saved out the 100 principal gradients feature matrix.
        % For other resolutions of PrincipalGrad, we need to load the whole matrix and get the first 59412*res rows.
        feature_mat_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'feature_mat', project_name, '100',...
            'input_corr_mat.mat'));
        idx = 59412*str2num(res);
        feature_names = fieldnames(feature_mat_path);
        params.feature_mat = feature_mat_path.(feature_names{2})(1:idx,:);
    else
        feature_mat_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'feature_mat', project_name, res,...
            'input_corr_mat.mat'));
        feature_names = fieldnames(feature_mat_path);
        params.feature_mat = feature_mat_path.(feature_names{2});
    end

    % Set up covariates matrix
    covariates_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params', 'regressors.mat'));
    params.covariates = covariates_path.regressors;

    % Set up target matrix
    y_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params', 'y_58_behaviors_3_components.mat'));
    params.y = y_path.y_in;

    % Set up output directory
    params.outdir = fullfile(out_dir, project_name, res, nsplit);

    % Set up output file name
    params.outstem = '58_behaviors_3_components';

    % Set up parameters for KRR: number of innerloop cross-validation folds
    params.num_inner_folds = 20;

    % Set up parameters for KRR: range of grid search for hyperparameter
    lambda_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params', 'lambda_set.mat'));
    params.lambda_set = lambda_path.lambda_set;

    ker_param_path = matfile(fullfile(REPDATA_DIR, 'HCP_data', 'input_params', 'ker_param.mat'));
    params.ker_param = ker_param_path.ker_param;

    params.metric = 'predictive_COD';

    params.with_bias = 1;

    % Call KRR workflow
    CBIG_KRR_workflow_LITE(params);
else
    disp('File already generated! Skip KRR workflow.')
end

end



