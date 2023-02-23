function CBIG_compute_multiKRR_PFM(feature_file, multiKRR_dir, sub_fold_file, score_ind, outstem, outdir)

% CBIG_compute_multiKRR_PFM(feature_file, multiKRR_dir, sub_fold_file, score_ind, outstem, outdir)
%
% This function computes the predictive network features of multi-kernel
% ridge regression (under <CBIG_CODE_DIR>/utilities/matlab/predictive_models/MultiKRR/) 
% by computing the covariance between the features and the predicted behavior (Haufe et al. 2014).
%
% Inputs:
%   - feature_file
%     The feature file used in the kernel regression (dim: #features x
%     #subjects). The variable storing the #features x #subjects matrix 
%     can be any name.
%
%   - multiKRR_dir
%     The directory where the multi-kernel ridge regression results are 
%     stored. This is the same directory you used as input argument 
%     `outdir` when you perform the multi-kernel ridge regression
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
%     from 1 to # Target Variables in your multi-kernel ridge regression
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
load(feature_file);
num_kernel = length(feature_files);
FC_all = CBIG_PFM_load_mat(feature_files{1});
[N_edge,N_sub] = size(FC_all);

features_all = zeros(N_edge,N_sub,num_kernel);
for i = 1:num_kernel
    FC_all = CBIG_PFM_load_mat(feature_files{i});
    FC_all = normalize(FC_all,'center');
    FC_all = normalize(FC_all,'norm',2);
    features_all(:,:,i) = FC_all;
    clear FC_all;
end

%% compute feature similarity matrix (FSM) 
FSM_all = zeros(N_sub,N_sub,num_kernel);
for i = 1:num_kernel
   
    tmp_features_all = features_all(:,:,i);
    
    tmp_FSM_all = CBIG_self_corr(tmp_features_all);
    if(~isempty(isnan(tmp_features_all)))
        nan_set = find(isnan(sum(tmp_features_all,1)));
        for nan_s = nan_set
            for all_s = 1:size(tmp_features_all,2)
                nan_mask = ~isnan(tmp_features_all(:,nan_s)) & ~isnan(tmp_features_all(:,all_s));
                tmp_FSM_all(nan_s,all_s) = CBIG_corr(tmp_features_all(nan_mask,nan_s),tmp_features_all(nan_mask,all_s));
                tmp_FSM_all(all_s,nan_s) = tmp_FSM_all(nan_s,all_s);
            end
        end
    end
    
    FSM_all(:,:,i) = tmp_FSM_all;
    clear tmp_features_all tmp_FSM_all;
end


%% compute activation
PFM_all_folds = zeros(N_edge*num_kernel,N_fold);
for i = 1:N_fold
    % load data
    load([multiKRR_dir '/optimal_acc/fold_' num2str(i) '/acc' outstem]);
    load([multiKRR_dir '/y/fold_' num2str(i) '/y_regress' outstem]);
    lambda_vect = opt_hyp{score_ind};
    test_ind = sub_fold(i).fold_index;
    train_ind = ~test_ind;
    
    y_train_resid = y_resid(train_ind,score_ind);
    FSM_train = FSM_all(train_ind,train_ind,:);
    FSM_pred = FSM_all(test_ind,train_ind,:);
        
    % training
    num_train_subj = size(FSM_train,1);
    y_tmp = repmat(y_train_resid,num_kernel,1);
    K_tmp = repmat(reshape(FSM_train,num_train_subj,num_train_subj*num_kernel),...
        num_kernel,1);
    kr_eye_tmp = eye(num_kernel*num_train_subj);
    kr_lambda_tmp = reshape(repmat(lambda_vect,num_train_subj,1),...
        num_train_subj*num_kernel,1);
    kr_tmp = bsxfun(@times, kr_eye_tmp, kr_lambda_tmp);
    K_tmp = K_tmp + kr_tmp;
    
    X = ones(num_train_subj*num_kernel,1);
    beta = (X' * (K_tmp\ X)) \ X' * (K_tmp \ y_tmp);
    alpha = K_tmp \ (y_tmp - X * beta);
    
    % prediction
    y_predicted_train = reshape(FSM_train,num_train_subj,...
        num_train_subj*num_kernel)*alpha + ones(num_train_subj,1) .* beta;
    
    % compute predictive-feature matrices
    features_train = zeros(N_edge*num_kernel,num_train_subj);
    for curr_kernel = 1:num_kernel
        features_train(1+(curr_kernel-1)*N_edge:curr_kernel*N_edge,:)...
            = features_all(:,train_ind,curr_kernel);
    end
    PFM_all_folds(:,i) = CBIG_PFM_cov_matrix(features_train',y_predicted_train)/std(y_predicted_train);
end

save([outdir '/PFM_score' num2str(score_ind) '_all_folds.mat'],'PFM_all_folds');
end


