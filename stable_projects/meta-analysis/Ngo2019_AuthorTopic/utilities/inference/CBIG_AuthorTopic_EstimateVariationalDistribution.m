function [newPhi, newVarBound] = ...
  CBIG_AuthorTopic_EstimateVariationalDistribution(params, ...
  newVarBound, varBoundFlag)
% [newPhi, newVarBound] = ...
%   CBIG_AuthorTopic_EstimateVariationalDistribution(params, ...
%   newVarBound, varBoundFlag)
%
% Update the parameters of the variation distribution phi of the
% Collapsed Variational Bayes (CVB) algorithm.
% Optionally update the lower bound of the data log likelihood.
%
% Input:
%  - params      : struct containing parameters of the author-topic model
%                  and the CVB algorithm.
%  - newVarBound : current estimate of the lower bound of the data log
%                  likelihood.
%  - varBoundFlag: if true, update the lower bound and do nothing otherwise.
% Output:
%  - newPhi      : updated estimate of the parameters of the variational
%                  distribution.
%  - newVarBound : updated estimate of the variational log likelihood.
%
% Example:
%   [newPhi, newVarBound] = ...
%     CBIG_AuthorTopic_EstimateVariationalDistribution(params, ...
%     newVarBound, true)
%   Update the parameters of the variational distribution and the lower
%   bound of the data log likelihood.
%
% Reference:
%   1) Eq(14) of the Supplemental S1 of Ngo et al., 2019,
%      "Beyond Consensus: Embracing Heterogeneity in Curated
%       Neuroimaging Meta-Analysis"
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  newPhi = cell(params.E, 1);
  for e = 1:params.E
    Ne = params.fociCounts(e);
    Te = params.expByNumTasks(e);
    newPhi{e} = zeros(Ne, params.K, Te, 'single');
  end;

  [newPhi, newVarBound] = CBIG_AuthorTopic_ComputeVariationalTerm__N_c(params, newPhi, newVarBound, varBoundFlag);
  [newPhi, newVarBound] = CBIG_AuthorTopic_ComputeVariationalTerm__N_cv(params, newPhi, newVarBound, varBoundFlag);
  [newPhi, newVarBound] = CBIG_AuthorTopic_ComputeVariationalTerm__N_t(params, newPhi, newVarBound, varBoundFlag);
  [newPhi, newVarBound] = CBIG_AuthorTopic_ComputeVariationalTerm__N_tc(params, newPhi, newVarBound, varBoundFlag);

  % normalize new phi
  newPhi = CBIG_AuthorTopic_NormalizePhi(newPhi, params);
