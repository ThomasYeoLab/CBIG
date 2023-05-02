function CBIG_ICCW_compute_krr_weights(krr_input,krr_result_dir,score_ind,outdir)

% function CBIG_ICCW_compute_krr_weights(krr_input,krr_result_dir,score_ind,outdir)
% 
% This function computes original regression weights for KRR. The original weights
% are a #features x #folds matrix which is saved into a mat file. 
% 
% Inputs:
%   - krr_input
%     Full path of the feature file used to calculate the kernels. A matrix
%     "feature_mat" is assumed to be saved in this file. This is normally the 
%     FC matrix that was used for prediction
%
%   - krr_result_dir
%     This refers to the directory in which the KRR results were saved.
%  
%   - score_ind
%     Score index to calculate regression weights for. Score index is in the same
%     order in which y variables were extracted for KRR.
%  
%   - outdir
%     Output directory in which to save the regression weights.
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% input
if(isstr(score_ind))
    score_ind = str2num(score_ind);
end
N_fold = 252;

%% load features
load(krr_input);
[N_edge,~] = size(FC_all);

FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);
features_all = FC_all;

%% compute FSM
FSM_all = corr(features_all);

%% compute activation
weights_all_folds = zeros(N_edge,N_fold);
load(fullfile(krr_result_dir,'no_relative_5_fold_sub_list.mat'));
load(fullfile(krr_result_dir,'final_result_all_score'));
for i = 1:N_fold
    % load data
    load(fullfile(krr_result_dir, 'y', strcat('fold_', num2str(i)), 'y_regress_all_score.mat'));
    lambda = optimal_lambda(i,score_ind);
    train_ind = sub_fold(i).fold_index == 0;
    
    y_train_resid = y_resid(train_ind,score_ind);
    FSM_train = FSM_all(train_ind,train_ind);
    features_train = features_all(:,train_ind);
    % number of subjects
    num_train_subj = size(FSM_train,1);
    
    % training
    X = ones(num_train_subj,1);
    K_r = FSM_train+lambda*eye(num_train_subj);
    y = y_train_resid;
    beta = (X' * (K_r\ X)) \ X' * (K_r \ y);
    alpha = K_r \ (y - X * beta);
    % prediction
    weights = features_train*alpha;
    % compute activation
    weights_all_folds(:,i) = weights;
end

save(fullfile(outdir, strcat('weights_score', num2str(score_ind), '_all_folds.mat')),'weights_all_folds');
end
