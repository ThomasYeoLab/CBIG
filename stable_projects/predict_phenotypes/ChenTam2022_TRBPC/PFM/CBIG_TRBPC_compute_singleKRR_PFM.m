function CBIG_TRBPC_compute_singleKRR_PFM(feature_file, singleKRR_dir, sub_fold_file, score_ind, outdir)

% CBIG_TRBPC_compute_singleKRR_PFM(feature_file,singleKRR_dir,score_ind,outdir)
%
% This function computes the predictive network features of single-kernel
% ridge regression by compute the covariance between the features and the
% predicted behavior.
%
% Inputs:
%   - feature_file
%     The feature file used in the kernel regression
%
%   - singleKRR_dir
%     The directory where the single-kernel ridge regression results are 
%     stored. This is the same directory you used as input argument 
%     `outdir` when you perform the kernel regression
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
%     from 1 to # Target Variables in your kernel regression
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

%% load  and normalize features
FC_all = CBIG_TRBPC_load_mat(feature_file);
[N_edge,~] = size(FC_all);

FC_all = normalize(FC_all,'center');
FC_all = normalize(FC_all,'norm',2);
features_all = FC_all;

%% compute FSM
FSM_all = corr(features_all);

%% compute activation
PFM_all_folds = zeros(N_edge,N_fold);
load(fullfile(singleKRR_dir,'final_result_all_score'));
for i = 1:N_fold
    % load data
    load([singleKRR_dir '/y/fold_' num2str(i) '/y_regress_all_score.mat']);
    lambda = optimal_lambda(i,score_ind);
    test_ind = sub_fold(i).fold_index;
    train_ind = ~test_ind;
    
    y_train_resid = y_resid(train_ind,score_ind);
    FSM_train = FSM_all(train_ind,train_ind);
    FSM_pred = FSM_all(test_ind,train_ind);
        
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
