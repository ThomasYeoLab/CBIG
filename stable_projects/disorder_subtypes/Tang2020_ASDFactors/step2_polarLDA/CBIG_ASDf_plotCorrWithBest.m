function CBIG_ASDf_plotCorrWithBest(inputDir, k, r_final, corr, outputDir)
% CBIG_ASDf_plotCorrWithBest(inputDir, k, r_final, corr, outputDir)
%
% This function plots the correlation between final estimate and other
% solutions, sorted from highest to lowest likelihood.
% 
% Input:
%     - inputDir:
%           Absolute path to the input directory where the model estimation
%           results are saved
%     - k:
%           A string indicating the number of factors
%     - r_final:
%           The final chosen estimate. E.g., in our paper, the solution
%           having the highest average correlation with all other solutions
%           was chosen as the final estimate.
%     - corr:
%           Rx1 cell array, where R is the number of random intializations.
%           Each cell is a kxR vector, where k is the number of factors.
%           corr{j,1} is the correlations between the j-th run and all
%           other runs. You can obtain correlation between all runs using
%           the function CBIG_ASDf_computCorrBetweenRuns.m
%     - outputDir (optional):
%           Absolute path to the output directory to save the plot
% 
% Example:
%	CBIG_ASDf_plotCorrWithBest('~/example_output/estimate','3','94',corr_k3,'~/example_output/visualizeFactors')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


%Check input variables
for i = 1: length(corr)
    tmp = corr{i};
    if size(tmp,1) ~= str2double(k) || size(tmp,2) ~= length(corr)
        error('Each cell in ''corr'' should have a dimension kxR');
    end
end

nRuns = length(corr);

% Rank the runs from highest to lowest likelihood
logLike = zeros(nRuns, 1);
for idx = 1:nRuns
    logLikes = load(fullfile(inputDir,sprintf('k%s/r%s', k, num2str(idx)),'likelihood.dat'));
    logLike(idx) = logLikes(end, 1);
end
[~, ind_h2l] = sort(logLike, 'descend');

ind_r_final = find(ind_h2l==str2double(r_final));

avgCorr = mean(corr{str2double(r_final)}, 1); % Get the average correlation between final estimate and other solutions
avgCorr = avgCorr(ind_h2l); % Sort correlations according to likelihood

% Plot
figure;
plot(1:nRuns, avgCorr, '-o');
hold on;
plot(ind_r_final, avgCorr(ind_r_final), '-o', 'MarkerEdgeColor','r');
box off;
xlabel('Random Initializations (highest likelihood to lowest)');
ylabel('Correlation with the Best Run');
xlim([1, nRuns]);
ylim([0, 1]);
set(gca, 'YTick', 0:0.05:1);
set(gca, 'TickDir', 'out');
grid on;

% Save the plot
if nargin > 4 && ~isempty(outputDir)
    output_name = fullfile(outputDir, ['k' k], ['r' r_final], 'corrWithAllRuns');
    hgexport(gcf, output_name);
    eps2xxx([output_name '.eps'], {'png'});
end
