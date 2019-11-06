function [ feat_train, feat_test ] = CBIG_FC_FeatSel( FC_train, FC_test, y_train, threshold )

% [ feat_train, feat_test ] = CBIG_FC_FeatSel( FC_train, FC_test, y_vec, fold_index, threshold )
% 
% This function selects a subset of functional connectivity (FC) entries as
% features for future predictive models. The selection is based on the
% correlation between FC and the target measure across the training
% subjects.
% 
% Inputs:
%   - FC_train: 
%     A 2D or 3D matrix containing functional connectivity of all training
%     subjects. If FC_train is 2D, it is assumed that FC_train is a K by N
%     matrix where K is the number of features and N is the number of 
%     training subjects. If FC_train is 3D, it is assumed to be a PxPxN
%     matrix, where P refers to the number of ROIs and N refers to the
%     number of training subjects. The 3D matrix will be converted to a 2D
%     matrix by taking the lower triangle of the 3D matrix.
% 
%   - FC_test:
%     A 2D or 3D matrix containing functional connectivity of all test
%     subjects. The dimension assumption is the same as FC_train.
% 
%   - y_train:
%     A column vector of the target measure to predict (e.g. fluid
%     intelligence). The length of y_train is the number of training
%     subjects.
% 
%   - threshold:
%     A scalar between 0 to 1. It is the percentage of total entries that
%     can be selected as featres.
% 
% Outputs:
%   - feat_train:
%     A 2D matrix of selected training features (dim: #features x #training
%     subjects). Every feature is standardized.
% 
%   - feat_test:
%     A 2D matrix of selected test features (dim: #features x #test
%     subjects). Every feature is standardized by using the mean and
%     standard deviation from the training set.
% 
% Written by Jingwei Li and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if(ndims(FC_train) == 3)
    temp = ones(size(FC_train,1), size(FC_train,2));
    tril_ind = tril(temp, -1);
    
    FC_train = reshape(FC_train, size(FC_train,1)*size(FC_train,2), size(FC_train, 3));
    FC_train = FC_train(tril_ind==1, :);
end
if(ndims(FC_test) == 3)
    temp = ones(size(FC_test,1), size(FC_test,2));
    tril_ind = tril(temp, -1);
    
    FC_test = reshape(FC_test, size(FC_test,1)*size(FC_test,2), size(FC_test,3));
    FC_test = FC_test(tril_ind==1, :);
end

train_corr = CBIG_corr(FC_train', y_train);

[~, sort_ind] = sort(abs(train_corr(:)), 'descend');
feature_ind = sort_ind(1:round(length(train_corr(:))*threshold));

feat_train = FC_train(feature_ind, :);
feat_test = FC_test(feature_ind, :);

%% normalization
mean_train = mean(feat_train, 2);
std_train = std(feat_train, 1, 2);

feat_train = bsxfun(@minus, feat_train, mean_train);
feat_train = bsxfun(@times, feat_train, 1./std_train);

feat_test = bsxfun(@minus, feat_test, mean_train);
feat_test = bsxfun(@times, feat_test, 1./std_train);

end

