function CBIG_compute_singleKRR_PFM(feature_file, singleKRR_dir, sub_fold_file, score_ind, outstem, outdir)

% CBIG_compute_singleKRR_PFM(feature_file, singleKRR_dir, sub_fold_file, score_ind, outstem, outdir)
%
% This function computes the predictive network features of single-kernel
% ridge regression (under <CBIG_CODE_DIR>/utilities/matlab/predictive_models/KernelRidgeRegression/) 
% by computing the covariance between the features and the predicted behavior (Haufe et al. 2014).
%
% Inputs:
%   - feature_file
%     The feature file used in the kernel regression (dim: #features x
%     #subjects). The variable storing the #features x #subjects matrix 
%     can be any name.
%
%   - singleKRR_dir
%     The directory where the single-kernel ridge regression results are 
%     stored. This is the same directory you used as input argument 
%     `outdir` when you perform the single-kernel ridge regression
%
%   - sub_fold_file
%     Full path of the cross-validation data split file. A structure
%     "sub_fold" is assumed to be saved in this file. 
%     sub_fold(i).fold_index is a #subjects x 1 binary vector. 1 refers to
%     the corresponding subject is a test subject in the i-th test fold. 0
%     refers to the corresponding subject is a training subject in the i-th
%     test fold.
%
%   - score_ind
%     A scalar. The index of the score you want to compute the PFM. Range
%     from 1 to # Target Variables in your single-kernel ridge regression
%
%   - outstem
%     A string that was used to describe the .mat file of the target
%     variable y after regression.
%
%   - outdir
%     Output directory for the predictive-feature matrices
%
% Outputs:
%   One # features by # folds predictive-feature matrix will be saved to the 
%   output directory for each behavior score
%
% Written by J Chen and N Wulan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if ~exist(outdir,'dir')
    mkdir(outdir)
end

if ~isempty(outstem)
    outstem = ['_' outstem];
end

%% get the cross-validation splits
if(isstr(score_ind))
    score_ind = str2num(score_ind);
end
load(sub_fold_file);
N_fold = length(sub_fold);

%% load and normalize features
FC_all = CBIG_PFM_load_mat(feature_file);
[N_edge,~] = size(FC_all);

FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);
features_all = FC_all;

%% compute feature similarity matrix (FSM)
FSM_all = CBIG_self_corr(features_all);
if(~isempty(isnan(features_all)))
    nan_set = find(isnan(sum(features_all,1)));
    for nan_s = nan_set
        for all_s = 1:size(features_all,2)
            nan_mask = ~isnan(features_all(:,nan_s)) & ~isnan(features_all(:,all_s));
            FSM_all(nan_s,all_s) = CBIG_corr(features_all(nan_mask,nan_s),features_all(nan_mask,all_s));
            FSM_all(all_s,nan_s) = FSM_all(nan_s,all_s);
        end
    end
end

%% compute activation
PFM_all_folds = zeros(N_edge,N_fold);
load(fullfile(singleKRR_dir, ['final_result', outstem]));
for i = 1:N_fold
    % load data
    load([singleKRR_dir '/y/fold_' num2str(i) '/y_regress' outstem]);
    lambda = optimal_lambda(i,score_ind);
    test_ind = sub_fold(i).fold_index;
    train_ind = ~test_ind;
    
    y_train_resid = y_resid(train_ind,score_ind);
    FSM_train = FSM_all(train_ind,train_ind);
        
    % training
    num_train_subj = size(FSM_train,1);
    X = ones(num_train_subj,1);
    K_r = FSM_train+lambda*eye(num_train_subj);
    y = y_train_resid;
    beta = (X' * (K_r\ X)) \ X' * (K_r \ y);
    alpha = K_r \ (y - X * beta);
    
    % prediction and compute predictive feature matrices
    y_predicted_train = FSM_train*alpha + ones(num_train_subj,1) .* beta;
    features_train = features_all(:,train_ind);
    PFM_all_folds(:,i) = CBIG_PFM_cov_matrix(features_train',y_predicted_train)/std(y_predicted_train);
end

save([outdir '/PFM_score' num2str(score_ind) '_all_folds.mat'],'PFM_all_folds');
end

