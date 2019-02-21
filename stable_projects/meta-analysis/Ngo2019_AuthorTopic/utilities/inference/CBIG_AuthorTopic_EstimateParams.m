function [theta, beta] = CBIG_AuthorTopic_EstimateParams(params)
% [theta, beta] = CBIG_AuthorTopic_EstimateParams(params)
%
% Compute estimates of the author-topic model's parameters (theta and
% beta) using the current estimate of the parameters of the variational
% distribution produced by the Collapsed Variational Bayes (CVB)
% algorithm.
%
% Input:
%  - params: struct containing parameters of the author-topic model and
%            variables estimated by the CVB algorithm.
% Output:
%  - theta : probability of a task recruiting a component
%            (Pr(component | task)).
%  - beta  : probability of a component activating a voxel
%            (Pr(voxel | component)).
%
% Example:
%   [theta, beta] = CBIG_AuthorTopic_EstimateParams(params)
%   Update the estimates of the theta, the probability of a task
%   recruiting a component (Pr(component | task)), and beta, probability
%   of a task recruiting a component (Pr(component | task)) using the
%   parameters and variables estimated by the CVB algorithm saved in
%   struct "params".
%
% Reference:
%   1) Eq(16) and (17) of the Supplemental S1 of Ngo et al., 2019,
%   "Beyond Consensus: Embracing Heterogeneity in Curated
%    Neuroimaging Meta-Analysis"
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  theta   = ones(params.T, params.K, 'single') * params.alpha;
  beta    = ones(params.K, params.V, 'single') * params.eta;

  for e = 1:params.E
    fociIndices    = params.fociIndices{e};
    taskIndices    = params.expByTask(e,:);
    Ad = params.expByNumTasks(e);

    tmp_1 = squeeze(sum(params.phi{e}, 1));
    if Ad > 1
        tmp_1 = tmp_1';
    end;
    theta(taskIndices, :) = theta(taskIndices, :) + tmp_1;

    tmp_2 = squeeze(sum(params.phi{e}, 3));
    beta(:, fociIndices) = beta(:, fociIndices) + tmp_2';
  end;

  theta = bsxfun(@times, theta, 1./sum(theta, 2));
  beta = bsxfun(@times, beta, 1./sum(beta, 2));
