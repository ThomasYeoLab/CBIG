function [corr, avg_corr_sorted, I] = CBIG_ASDf_computeCorrBetweenRuns(inputDir, k)
% [corr, avg_corr_sorted, I] = CBIG_ASDf_computeCorrBetweenRuns(inputDir, k)
% 
% This function computes pair-wise correlation between each run (random
% initialization). In addition, for reach run, it computes the average correlation
% with all other runs, and sort the runs according to the average
% correlation.
%
% Input:
%     - inputDir:
%           Absolute path to the directory where the results of all the random initializations are saved
%     - k:
%           A string indicating the number of factors, e.g., k = '3'
% Output:
%     - corr:
%           A Nx1 cell array, where the i-th row is a kxN matrix (N:number of
%           runs; k: number of factors) representing the correlation between
%           the i-th run and all the runs (factors matched using Hungarian
%           matching method).
%     - avg_corr_sorted:
%           A 1xN array, where the i-th entry is the average correlation between
%           the i-th run and all the runs, sorted from largest to smallest.
%     - I:
%           1xN array. Run number sorted from largest average correlation to smallest
%           average correlation (i.e., run number corresponding to
%           avg_corr_sorted).
%
% Example:
%	[corr, avg_corr_sorted, I] = CBIG_ASDf_computeCorrBetweenRuns('~/example_output/estimate','3');
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

dirList = dir([inputDir sprintf('/k%s/r*', k)]);
nRuns = numel(dirList);

corr = cell(nRuns,1);

for i = 1:nRuns
    for j = (i+1):nRuns
%         fprintf('Correlation between run%d and run%d:\n', i, j);
        corr{i}(:,j) = CBIG_ASDf_corrTwoRuns(inputDir, inputDir, k, num2str(i), num2str(j));
    end
end


for i = 1:nRuns
    corr{i}(:,i) = ones(str2double(k),1);
    for j = 1:(i-1)
        corr{i}(:,j) = corr{j}(:,i);
    end
end

for i = 1:length(corr)
    avg_corr(i) = mean(mean(corr{i},2));
end

[avg_corr_sorted, I] = sort(avg_corr,'descend');
% save([outputDir '/corrBetweenRuns_k=' k], 'corr', 'avg_corr_sorted', 'I');
