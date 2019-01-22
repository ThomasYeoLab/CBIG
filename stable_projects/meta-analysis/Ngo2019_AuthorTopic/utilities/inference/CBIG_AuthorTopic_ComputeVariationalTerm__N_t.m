function [newPhi, newVarBound] = ...
  CBIG_AuthorTopic_ComputeVariationalTerm__N_t(params, newPhi, ...
  newVarBound, varBoundFlag)
% [newPhi, newVarBound] = ...
% CBIG_AuthorTopic_ComputeVariationalTerm__N_t(params, newPhi, ...
% newVarBound, varBoundFlag)
%
% Compute the expectation and variance of the the count N_{.t..}^{-ef}
% in the update equation for the parameter of the variation distribution
% phi.
% Optionally update the variational bound of the data log likelihood.
%
% Input:
%  - params      : struct containing parameters of the author-topic model
%                  and the Collapsed Variational Bayes (CVB) algorithm.
%  - newPhi      : current estimate of the paramters of the variational
%                  distribution.
%  - newVarBound : current estimate of the variational bound of the data
%                  log likelihood.
%  - varBoundFlag: if true, update the lower bound and do nothing
%                  otherwise.
% Output:
%  - newPhi      : updated estimate of the parameters of the variational
%                  distribution.
%  - newVarBound : updated estimate of the variational bound.
%
% Example:
%   [newPhi, newVarBound] = ...
%     CBIG_AuthorTopic_ComputeVariationalTerm__N_t(params, newPhi, ...
%     newVarBound, true)
% Update the parameters of the variational distribution and the
% variational bound of the data log likelihood by computing the
% expectation and variance of the count N_{.t..}^{-ef}.
%
% Reference:
%   1) Eq(14) of the Supplemental S1 of Ngo et al., 2019,
%   "Beyond Consensus: Embracing Heterogeneity in Neuroimaging
%    Meta-Analysis"
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_3                           = cell(params.E, 1);
  Var_3                         = cell(params.E, 1);
  sum_E_3                       = zeros(1, params.T, 'single');
  sum_Var_3                     = zeros(1, params.T, 'single');
  phi_E_3                       = cell(params.E, 1);
  phi_Var_3                     = cell(params.E, 1);

  for e = 1:params.E
    taskIndices                = params.expByTask(e,:);
    % N_{.t..}
    % for E_q(N_{.t..})
    phi_E_3{e}                  = squeeze(sum(params.phi{e},2));
    e_tmp_1                     = bsxfun(@times, phi_E_3{e}, params.Fc{e});
    e_tmp_2                     = squeeze(sum(e_tmp_1, 1));
    sum_E_3(taskIndices)     = sum_E_3(taskIndices) + e_tmp_2;
    clear e_tmp_1 e_tmp_2;
    % for Var_q(N_{kw_{en}})
    phi_Var_3{e}                = phi_E_3{e} .* (1 - phi_E_3{e});
    var_tmp_1                   = bsxfun(@times, phi_Var_3{e}, params.Fc{e});
    var_tmp_2                   = squeeze(sum(var_tmp_1, 1));
    sum_Var_3(taskIndices)   = sum_Var_3(taskIndices) + var_tmp_2;
    % clear
    clear var_tmp_1 var_tmp_2;
  end;

  if varBoundFlag
    % Compute the variational bound
    for t = 1:params.T
      n_3               = single([0:(uint32(params.Nt(t)) - 1)]);
      log_tmp_3         = log(single(params.K * params.alpha + [1:params.Nt(t)] - 1));
      mu_3              = single(sum_E_3(t));
      sigma_3           = single(sqrt(sum_Var_3(t)));

      if sigma_3 <= 0
        newVarBound   = newVarBound - sum(log_tmp_3(1:floor(mu_3)));
      else
        q_3             = 1 - normcdf(n_3, mu_3, sigma_3);
        newVarBound   = newVarBound - sum(q_3 .* log_tmp_3);
      end;

      clear mu_3 sigma_3;
    end;
    clear n_3 log_tmp_3
  end;


  % Deduct variables of d, n
  for e = 1:params.E
    Ne                          = params.fociCounts(e);
    taskIndices                = squeeze(params.expByTask(e,:));
    % N_{.t..}
    % for E_q(N_{.t..})
    tmp_1                       = repmat(sum_E_3(taskIndices), [Ne 1]);
    E_3{e}                      = tmp_1 - phi_E_3{e};
    clear tmp_1;
    % for Var_q(N_{.t..})
    tmp_2                       = repmat(sum_Var_3(taskIndices), [Ne 1]);
    Var_3{e}                    = tmp_2 - phi_Var_3{e};
    % clear
    clear tmp_2
  end;
  clear sum_E_3 sum_Var_3 phi_E_3 phi_Var_3

  % update phi
  for e = 1:params.E
    X_3_tmp                     = params.K * params.alpha + E_3{e};
    X_3                         = log(X_3_tmp);
    X_3                         = permute(X_3, [1 3 2]);
    Y_3_tmp                     = Var_3{e} ./ (X_3_tmp.^2);
    Y_3                         = permute(Y_3_tmp, [1 3 2]);

    newPhi{e}                  = bsxfun(@minus, newPhi{e}, X_3);
    newPhi{e}                  = bsxfun(@plus, newPhi{e}, 0.5 * Y_3);

    clear X_3 X_3_tmp Y_3 Y_3_tmp
  end;
  clear E_3 Var_3
