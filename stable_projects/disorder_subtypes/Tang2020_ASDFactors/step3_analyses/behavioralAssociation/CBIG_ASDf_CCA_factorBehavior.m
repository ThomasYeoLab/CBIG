function [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior(ind_scores,...
    id_keep, reg, PAPset, path_subInfo, path_behav_scores, path_factorLoading, factor_order,...
    factor_idx)
% [A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] =
% CBIG_ASDf_CCA_factorBehavior(ind_scores, id_keep, reg, PAPset, path_subInfo, 
% path_behav_scores, path_factorLoading, factor_order, factor_idx)
% 
% This function performs CCA between factor loading and behavioral scores
% with permutation test (within site permutation).
%
% Input:
%      - ind_scores:
%           1xQ vector, where Q is the number of behavioral scores. Column
%           indices of behavioral scores in the .csv file specified by
%           path_behav_scores
%      - id_keep:
%           Nx1 cell array. IDs of subjects of interest
%      - reg:
%           NxV matrix, where V is the number of regressors
%      - PAPset:
%           N x Nperm matrix, where Nperm is the number of permutations. 
%           Permutation set generated using PALM package. 
%           The 1st column should be the original order (non-permuted).
%      - path_subInfo:
%           .csv file containing all subjects' demographic information, e.g., diagnosis, age, sex etc.
%      - path_behav_scores:
%           .csv file containing all behavioral scores
%      - path_factorLoading:
%           .txt file containing all ASD subjects' factor compositions
%      - factor_order:
%           Order of the factors. E.g., [1 3 2] means that the factors will 
%           be re-ordered so that the 1st factor indicates factor 1, 3rd
%           factor indicates factor 2, and 2nd factor indicates factor 3.
%      - factor_idx:
%           Index of the factor
%
% Output:
%      - A:
%           Qx1 matrix, where Q is the number of behavioral scores. This is
%           the canonical coefficients for behavioral scores
%      - strucCorr_score:
%           Qx1 matrix, where Q is the number of behavioral scores. This is 
%           the structure coefficients for behavioral scores. i.e., the
%           Pearson's correlation between canonical coefficients A and the
%           raw behavioral scores.
%      - B:
%           The canonical coefficient for factor loading.
%      - U:
%           Mx1 matrix, where M is the number of ASD participants. This is 
%           the canonical covariate (mode) for behavioral scores.
%      - V:
%           Mx1 matrix, where M is the number of ASD participants. This is the
%           canonical covariate (mode) for factor loading.
%      - R:
%           The canonical correlation (Pearson's correlation between U and V).
%           This is what CCA maximizes.
%      - Rp:
%           Nperm x 1 matrix, where Nperm is the number of permutations. 
%           This is the canonical correlations of all permutations.
%      - pVal:
%           p-value from permutation test.
%      - Ncca:
%           The number of significant CCA mode survived after permutation test.
%
%
% Example:
%	[A, strucCorr_score, B, U, V, R, Rp, pVal, Ncca] = CBIG_ASDf_CCA_factorBehavior([1,2,3], 
%   id_keep, reg, PAPset, '../examples/input/subInfo_est.csv', '../examples/input/behavior_scores_est.csv', 
%   '~/example_output/visualizeFactors/k3/r94/factorComp.txt', [3 2 1], 1)
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Check input variable
if  size(ind_scores,1) > 1
    error('Input argument ''ind_scores'' should be a row vector');
end

if size(factor_order,1) > 1
    error('Input argument ''factor_order'' should be a row vector');
end

%% Get behavioral scores & factor loadings
[~, ~, ~, ~, ~, ~, id_scores, id_factorLoading] = CBIG_ASDf_getSubData(path_subInfo, ...
    id_keep, path_behav_scores, ind_scores, path_factorLoading);

scores = cell2mat(id_scores(:,2:end));
factorLoading = cell2mat(id_factorLoading(:,2:end));

factorLoading = factorLoading(:,factor_order); % Reorder the factors
fl = factorLoading(:,factor_idx);

conf = reg;

Nperm = size(PAPset,2);
%% Regress out age, sex, motion and sites from factor loading
fl_d = fl-conf*(pinv(conf)*fl); % Regress out regressors from factor loading

%% Regress out age, sex, motion and sites from behavioral scores
scores_d = zeros(size(scores));
for i=1:size(scores,2)
    % Regress out regressors from behavioral scores
    scores_d(:,i) = scores(:,i)-conf*(pinv(conf)*scores(:,i)); 
end

%% CCA
% Age, sex, motion & sites regressed on both sides
[A, B, R, U, V, ~] = canoncorr(scores_d,fl_d); 

%% Permutation test
Rp=zeros(Nperm,size(factorLoading,2)); 
clear pVal;
for j=1:Nperm
    fprintf('Permutation No.%d\n',j);
    % Age, sex, motion & sites regressed on both sides
    [Ap,Bp,Rp(j,:),Up,Vp]=canoncorr(scores_d,fl_d(PAPset(:,j),:)); 
end

pVal = (1+sum(Rp(2:end,1)>=R))/Nperm;
pVal
Ncca = sum(pVal<0.05)

% Always keep B (canonical coefficient for factor loading) positive 
if B < 0
    A = -A;
    U = -U;
    B = -B;
    V = -V;
end

% Compute structure coefficient for behavioral scores
strucCorr_score(:,1) = corr(U(:,1),scores,'rows','complete'); 

