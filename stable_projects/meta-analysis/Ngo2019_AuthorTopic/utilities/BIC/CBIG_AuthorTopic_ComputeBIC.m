function CBIG_AuthorTopic_ComputeBIC(allKs, maskPath, bestSolutionDir, ...
    outputDir)
% CBIG_AuthorTopic_ComputeBIC(allKs, maskPath, bestSolutionDir, ...
%   outputDir)
%
% Compute Bayesian Information Criterion (BIC) measure of the author-topic
% model's parameter estimates across different numbers of components.
%
% Input:
%  - allKs          : vector containing the numbers of components.
%  - maskPath       : path to the mask defining the voxels that are used
%                     for BIC computation. This should be the mask of
%                     activation across all experiments included in the
%                     meta-analysis, produced by
%                     CBIG_AuthorTopic_GenerateCVBDataFromText.m
%  - bestSolutionDir: path to directory containing the best model estimates
%                     for each given number of components.
%  - outputDir      : directory containing all intermediate files and the
%                     final BIC measures.
% Output:
%  - The BIC measure of the author-topic model parameters estimate with a
%    given number of cognitive component K is stored at
%    <outputDir>/K<K>_BIC.mat. The plot of BIC measures across all numbers
%    of cognitive components is also produced.
%
% Example:
%   CBIG_AuthorTopic_ComputeBIC([1:5], '/Work/AT/expMask.nii.gz', ...
%     '/Work/AT/outputs/bestSolution', '/Work/AT/BIC')
%   Compute the BIC measures of the best estimates of the author-topic
%   model parameters for 1 to 5 cognitive components. The estimates are
%   from '/Work/AT/outputs/bestSolution/BestSolution_K001.mat' to
%   '/Work/AT/outputs/bestSolution/BestSolution_K005.mat'. The mask comes
%   from '/Work/AT/expMask.nii.gz'. The output BIC measures are saved under
%   '/Work/AT/BIC'.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  if ~exist(outputDir, 'dir')
    mkdir(outputDir);
  end

  componentSmoothnessDir = fullfile(outputDir, 'ComponentSmoothness');
  if ~exist(componentSmoothnessDir, 'dir')
    mkdir(componentSmoothnessDir);
  end

  allBICs = zeros(length(allKs), 1);
  allLogLikelihoods = zeros(length(allKs), 1);
  allComplexityCosts = zeros(length(allKs), 1);
  allFreeDimensionsCounts = zeros(length(allKs), 1);
  count = 1;

  CBIG_AuthorTopic_EstimateComponentSmoothness(allKs, maskPath, bestSolutionDir, componentSmoothnessDir);

  for K = allKs
    disp(['K: ' num2str(K)]);
    smoothness = load(fullfile(componentSmoothnessDir, ['K' num2str(K) '_smoothness.mat']));
    solution = load(fullfile(bestSolutionDir, ['BestSolution_K' num2str(K, '%03d') '.mat']));

    freeDimensionsCount = smoothness.totalReselsCount + size(solution.params.theta, 1) * (K-1);
    [bic, logLikelihood, complexityCost] = CBIG_AuthorTopic_ComputeBICFromATParams(solution.params, freeDimensionsCount);

    save(fullfile(outputDir, ['K' num2str(K) '_BIC.mat']), 'bic');
    allBICs(count) = bic;
    allLogLikelihoods(count) = logLikelihood;
    allFreeDimensionsCounts(count) = freeDimensionsCount;
    allComplexityCosts(count) = complexityCost;
    count = count + 1;
  end

  figure;
  plot(allKs, allBICs, '-ko', 'LineWidth', 2);
  title('Bayesian Information Criteria (BIC)', 'FontSize', 19);
  xlabel('# Components', 'FontSize', 19);
  ylabel('Bayesian Information Criterion (BIC)', 'FontSize', 19);
  set(gca, 'XTick', allKs);
  set(gca, 'FontSize', 19);
  set(gca, 'TickDir', 'out');
  box off;
  num_ticks = 4;
  L = get(gca, 'YLim');
  set(gca, 'YTick', linspace(L(1), L(2), num_ticks));
  saveas(gcf, fullfile(outputDir, 'BIC_3dFWHMx'), 'epsc');
