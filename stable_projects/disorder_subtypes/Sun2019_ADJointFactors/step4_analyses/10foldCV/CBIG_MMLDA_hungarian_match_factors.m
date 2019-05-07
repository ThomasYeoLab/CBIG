function [in_order, avg_corr] = CBIG_MMLDA_hungarian_match_factors(ref_dir, in_dir, ref_order)
% [in_order, avg_corr] = CBIG_MMLDA_hungarian_match_factors(ref_dir, in_dir, ref_order)
%
% Do a hungarian matching between reference betas and input betas and return the order and correlation of input betas.
%
% Input:
%   - ref_dir   : directory of reference "final.beta1", "final.beta2"
%   - in_dir    : directory of input "final.beta1", "final.beta2"
%   - ref_order : order of reference beta
%
% Output:
%   - in_order  : order of input beta
%   - avg_corr  : average correlation between input and reference beta
%
% Example:
%   [in_order, avg_corr] = CBIG_MMLDA_hungarian_match_factors(ref_dir, in_dir, [1 3 2])
%
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% load reference beta1 and beta2 and reorder it
ref_beta1 = exp(load([ref_dir '/final.beta1']));
ref_beta2 = exp(load([ref_dir '/final.beta2']));
ref_beta1 = ref_beta1(ref_order, :);
ref_beta2 = ref_beta2(ref_order, :);

% load input beta1 and beta2
in_beta1 = exp(load([in_dir '/final.beta1']));
in_beta2 = exp(load([in_dir '/final.beta2']));

%%% Reorder subtypes to obtain the maximal correlation coefficients
% Construct the COST matrix
%         pos1 pos2 ...
% topic1
% topic2
% ...
k = size(ref_beta1, 1);
costMat = zeros(k, k);
for rowIdx = 1:k
    for colIdx = 1:k
        % Assign beta (jobs, column) to beta_best (workers, row)
        corrMat = corrcoef(ref_beta1(rowIdx, :)', in_beta1(colIdx, :)');
        costMat(rowIdx, colIdx) = 1-corrMat(1, 2);
    end
end

% Run the Hungarian matching algorithm
% order: each row (worker)'s matched column (job)
[in_order, ~] = munkres(costMat);

% Recompute the avergae correlation with sorted topics
corr1 = zeros(k, 1);
for idx = 1:k  
    corrMat = corrcoef(ref_beta1(idx, :)', in_beta1(in_order(idx), :)');
    corr1(idx) = corrMat(1, 2); 
end

corr2 = zeros(k, 1);
for idx = 1:k  
    corrMat = corrcoef(ref_beta2(idx, :)', in_beta2(in_order(idx), :)');
    corr2(idx) = corrMat(1, 2); 
end

avg_corr = mean([corr1; corr2]);