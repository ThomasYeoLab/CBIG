function order = CBIG_ASDf_hunMatch(K, Mean_ref, Mean)
% order = CBIG_ASDf_hunMatch(K, Mean_ref, Mean)
% 
% This function performs Hungarian matching between the factor-specific RSFC
% patterns (i.e., E(RSFC patterns|Factor)) of two random initializations.
%
% Input:
%     - K:
%           An integer indicating the number of latent factors
%     - Mean_ref:
%           E(RSFC patterns|Factor) of the reference solution
%     - Mean:
%           E(RSFC patterns|Factor) of the solution to be matched with
%           the reference solution
% Output:
%     - order:
%           Order of the latent factors that are matched with the reference
%           solution
% 
% Example:
%	order = CBIG_ASDf_hunMatch(3, mean_run1, mean_run10)
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

costMat = zeros(K, K);
for rowIdx = 1:K
    for colIdx = 1:K
        corrMat = corrcoef(Mean_ref(rowIdx, :)', Mean(colIdx, :)');
        costMat(rowIdx, colIdx) = 1-corrMat(1, 2);
    end
end

% Run the Hungarian matching algorithm
[order, ~] = munkres(costMat);
