function CBIG_compute_LRR_fitrlinear_PFM(feature_file, LRR_dir, sub_fold_file, score_ind, outstem, outdir)

% CBIG_compute_LRR_fitrlinear_PFM(feature_file, LRR_dir, sub_fold_file, score_ind, outstem, outdir)
%
% This function computes the predictive network features of linear
% ridge regression (under <CBIG_CODE_DIR>/utilities/matlab/predictive_models/LRR_fitrlinear/) 
% by computing the covariance between the features and the predicted behavior (Haufe et al. 2014).
%
% Inputs:
%   - feature_file
%     The feature file used in the kernel regression (dim: #features x
%     #subjects). The variable storing the #features x #subjects matrix 
%     can be any name.
%
%   - LRR_dir
%     The directory where the linear ridge regression results are 
%     stored. This is the same directory you used as input argument 
%     `outdir` when you perform the linear ridge regression
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
%     from 1 to # Target Variables in your linear ridge regression
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

normalized_FC_all = normalize(FC_all,'center');
normalized_FC_all = normalize(normalized_FC_all,'norm',2);

%% compute activation
PFM_all_folds = zeros(N_edge,N_fold);
for i = 1:N_fold
    % load data
    load([LRR_dir '/y/fold_' num2str(i) '/y_regress' outstem]);
    load(fullfile(LRR_dir, 'params', ['fold_' num2str(i)], ['selected_parameters' outstem]));
    lambda = curr_lambda(score_ind);
    threshold = curr_threshold(score_ind);
    test_ind = sub_fold(i).fold_index;
    train_ind = ~test_ind;
    
    y_train_resid = y_resid(train_ind,score_ind);
    features_train = FC_all(:,train_ind);
    feature_test = FC_all(:,test_ind);
    [features_train, ~] = CBIG_FC_FeatSel( features_train, feature_test, y_train_resid, threshold);
    normalized_features_train = normalized_FC_all(:,train_ind);
        
    % training
    Mdl = fitrlinear(features_train,y_train_resid,'ObservationsIn','columns', 'Lambda',lambda, 'Learner',...
    'leastsquares', 'Regularization','ridge');

    % prediction and compute predictive feature matrices
    y_predicted_train = predict(Mdl, features_train,'ObservationsIn','columns');
    PFM_all_folds(:,i) = CBIG_PFM_cov_matrix(normalized_features_train',y_predicted_train)/std(y_predicted_train);
end

save([outdir '/PFM_score' num2str(score_ind) '_all_folds.mat'],'PFM_all_folds');
end

