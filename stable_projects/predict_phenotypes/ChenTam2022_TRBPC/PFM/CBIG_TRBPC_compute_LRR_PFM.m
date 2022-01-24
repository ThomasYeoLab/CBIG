function CBIG_TRBPC_compute_LRR_PFM(feature_file, LRR_dir, sub_fold_file, score_ind, outdir)

% CBIG_TRBPC_compute_LRR_PNF(feature_file,LRR_dir,score_ind,outdir)
%
% This function computes the predictive network features of linear
% ridge regression by computing the covariance between the features
% and the predicted behavior.
%
% Inputs:
%   - feature_file
%     The feature file used in the linear regression
%
%   - LRR_dir
%     The directory where the LRR results are stored. This is the same
%     directory you used as input argument `outdir` when you perform the
%     linear ridge regression
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

load(fullfile(LRR_dir,'param.mat'));
N_fold = length(param.sub_fold);

feat_mat = param.FC_mat;
N_feature = size(feat_mat,1);
i = score_ind;
PFM_all_folds = zeros(N_feature,N_fold);
for j = 1:N_fold
    ind = param.sub_fold(j).fold_index;
    features_train = feat_mat(:,~ind);
    load(fullfile(LRR_dir,['behav_' num2str(i)],'params',['fold_' num2str(j)],'opt_results_all_score.mat'));
    PFM_all_folds(:,j) = cov_matrix(features_train',pred_train)/std(pred_train);
end
save([ outdir '/PFM_score' num2str(i) '_all_folds.mat'],'PFM_all_folds');
rmpath(paths)
end

function cov_xy = cov_matrix(x,y)
% x: N*K1 matrix, N is # observation, K is # variables
% y: N*K2 matrix, N is # observation, K is # variables
% cov_xy: K1*K2 covariance matrix
cov_xy = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);
end
