function [PFM] = CBIG_compute_PFM_general(feature_file, y_pred_train_file)

% [PFM] = CBIG_compute_PFM_general(feature_file, y_pred_train_file)
% 
% This function computes the predictive network features for general purpose
% by computing the covariance between the features and the predicted behavior
% (Haufe et al. 2014).
%
% Inputs:
%   - feature_file
%     The feature file used in predictive model (dim: #features x #subjects).
%     The variable storing the #features x #subjects matrix can be any name.
%
%   - y_pred_train_file
%     The predicted y for training subjects from the predictive model
%     (dim: #subjects x 1). The variable storing the #subjects x 1 vector 
%     can be any name.
%
% Outputs:
%   - PFM
%      The predictive-feature matrix (dim: #features x 1).
%
% Written by Naren Wulan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% load features
FC = CBIG_PFM_load_mat(feature_file);
FC = normalize(FC,'center');
FC = normalize(FC,'norm',2);

%% load predicted y
y_predicted_train = CBIG_PFM_load_mat(y_pred_train_file);

%% compute predictive feature matrices
PFM = CBIG_PFM_cov_matrix(FC',y_predicted_train)/std(y_predicted_train);
end

