function CBIG_TRBPC_compute_multiKRR_PFM( feature_file, multiKRR_dir, sub_fold_file, score_ind, outdir)

% CBIG_TRBPC_compute_multiKRR_PFM(feature_file,multiKRR_dir,score_ind,outdir)
%
% This function computes the predictive network features of multi-kernel
% ridge regression by compute the covariance between the features and the
% predicted behavior.
%
% Inputs:
%   - feature_file
%     The feature file used in the multi-kernel regression
%
%   - multiKRR_dir
%     The directory where the multiKRR results are stored. This is the same
%     directory you used as input argument `outdir` when you perform the
%     multi-kernel regression
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
%     A scalar. The index of the score you want to compute the PNF. Range
%     from 1 to # Target Variables in your multi-kernel regression
%
%   - outdir
%     Output directory for the predictive-feature matrices
%
% Outputs:
%   One # features by # folds predictive-feature matrix will be saved to the 
%   output directory for each behavior score
%
% Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

project_code_dir = fullfile(getenv('CBIG_CODE_DIR'),'stable_projects','predict_phenotypes', 'ChenTam2022_TRBPC');
paths = genpath(project_code_dir);
addpath(paths)

if ~exist(outdir,'dir')
    mkdir(outdir)
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
FC_all = CBIG_TRBPC_load_mat(feature_files{1});
[N_edge,N_sub] = size(FC_all);

features_all = zeros(N_edge,N_sub,num_kernel);
for i = 1:num_kernel
    FC_all = CBIG_TRBPC_load_mat(feature_files{i});
    FC_all = normalize(FC_all,'center');
    FC_all = normalize(FC_all,'norm',2);
    features_all(:,:,i) = FC_all;
    clear FC_all;
end

%% compute FSM
FSM_all = zeros(N_sub,N_sub,num_kernel);
for i = 1:num_kernel
    FSM_all(:,:,i) = corr(features_all(:,:,i));
end
%% compute activation
PFM_all_folds = zeros(N_edge*num_kernel,N_fold);
for i = 1:N_fold
    % load data
    load([multiKRR_dir '/optimal_acc/fold_' num2str(i) '/acc_all_score.mat']);
    load([multiKRR_dir '/y/fold_' num2str(i) '/y_regress_all_score.mat']);
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
    PFM_all_folds(:,i) = cov_matrix(features_train',y_predicted_train)/std(y_predicted_train);
end

save([outdir '/PFM_score' num2str(score_ind) '_all_folds.mat'],'PFM_all_folds');
rmpath(paths)
end

function cov_xy = cov_matrix(x,y)
% x: N*K1 matrix, N is # observation, K is # variables
% y: N*K2 matrix, N is # observation, K is # variables
% cov_xy: K1*K2 covariance matrix
cov_xy = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);
end
