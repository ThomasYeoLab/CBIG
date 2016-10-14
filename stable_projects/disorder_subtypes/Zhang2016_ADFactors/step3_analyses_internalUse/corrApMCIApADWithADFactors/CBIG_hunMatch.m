function avgCorr = CBIG_hunMatch(beta_best, beta)

% avgCorr = CBIG_hunMatch(beta_best, beta)
%%% Reorder subtypes to obtain the maximal correlation coefficients
% Construct the COST matrix
%         pos1 pos2 ...
% topic1
% topic2
% ...
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

K = size(beta, 1);
costMat = zeros(K, K);
for rowIdx = 1:K
    for colIdx = 1:K
        % Assign beta (jobs, column) to beta_best (workers, row)
        corrMat = corrcoef(beta_best(rowIdx, :)', beta(colIdx, :)');
        costMat(rowIdx, colIdx) = 1-corrMat(1, 2);
    end
end

% Run the Hungarian matching algorithm
% order: each row (worker)'s matched column (job)
[order, ~] = munkres(costMat);

% Recompute the avergae correlation with sorted topics
corr = zeros(K, 1);
for idx = 1:K  
	corrMat = corrcoef(beta_best(idx, :)', beta(order(idx), :)');
    corr(idx) = corrMat(1, 2); 
end
corr
avgCorr = mean(corr);