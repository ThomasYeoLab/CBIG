function [bic, logLikelihood, complexityCost] = ...
  CBIG_AuthorTopic_ComputeBICFromATParams(params, freeDimensionsCount)
% [bic, logLikelihood, complexityCost] = ...
%   CBIG_AuthorTopic_ComputeBICFromATParams(params, freeDimensionsCount)
%
% Compute Bayesian Information Criterion (BIC) measure of the estimates of
% the author-topic model's parameters produced by the Collapsed Variational
% Bayes (CVB) algorithm and the parameters' number of free dimensions.
%
% Input:
%  - params             : estimates of the author-topic model parameters
%                         produced by CVB algorithm.
%  - freeDimensionsCount: number of free dimensions of the parameter
%                         estimates.
% Output:
%  - bic           : BIC measure.
%  - logLikelihood : data log likelihood of the model parameter estimates.
%  - complexityCost: measure of the model parameter estimates' complexity.
%
% Example:
%   [bic, logLikelihood, complexityCost] = ...
%     CBIG_AuthorTopic_ComputeBICFromATParams(bestParams, ..
%     freeDimensionsNumber)
%   Compute BIC measure of the author-topic model's parameters saved in
%     bestParams and the parameters' number of free dimensions
%     freeDimensionsNumber.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  logLikelihood = 0;
  wordCountSum = 0;

  beta = params.beta;
  theta = params.theta;

  for e = 1:params.E
    wordIndices = logical(params.Fb(e, :));
    betaForDoc = beta(:, wordIndices);
    wordCount = params.Fc{e}';

    authorsForDoc = logical(full(params.expByTask(e, :)));
    authorsForDocCount = sum(authorsForDoc);
    thetaForDoc = theta(authorsForDoc, :);

    docLogLikelihood = betaForDoc' * thetaForDoc';
    docLogLikelihood = log(sum(docLogLikelihood, 2));
    docLogLikelihood = bsxfun(@times, wordCount', docLogLikelihood);
    docLogLikelihood = sum(docLogLikelihood) - log(authorsForDocCount) * sum(wordCount);

    logLikelihood = logLikelihood + docLogLikelihood;
    wordCountSum = wordCountSum + sum(wordCount);
  end

  bic = logLikelihood - 0.5 * freeDimensionsCount * log(wordCountSum);
  complexityCost = -0.5 * freeDimensionsCount * log(wordCountSum);
