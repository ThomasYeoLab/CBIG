function [newPhi, newVarBound] = ...
  CBIG_AuthorTopic_ComputeVariationalTerm__N_tc(params, newPhi, ...
  newVarBound, varBoundFlag)
% [newPhi, newVarBound] = ...
% CBIG_AuthorTopic_ComputeVariationalTerm__N_tc(params, newPhi, ...
% newVarBound, varBoundFlag)
%
% Compute the expectation and variance of the the count N_{.tc.}^{-ef}
% in the update equation for the parameter of the variation distribution
% phi.
% Optionally update the lower bound of the data log likelihood.
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
%  - newVarBound : updated estimate of te lower bound.
%
% Example:
%   [newPhi, newVarBound] = ...
%     CBIG_AuthorTopic_ComputeVariationalTerm__N_tc(params, newPhi, ...
%     newVarBound, true)
% Update the parameters of the variational distribution and the
% variational bound of the data log likelihood by computing the
% expectation and variance of the count N_{.tc.}^{-ef}.
%
% Reference:
%   1) Eq(14) of the Supplemental S1 of Ngo et al., 2019,
%      "Beyond Consensus: Embracing Heterogeneity in Neuroimaging
%       Meta-Analysis"
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_4                               = cell(params.E, 1);
  Var_4                             = cell(params.E, 1);
  sum_E_4                           = zeros(params.K, params.T, 'single');
  sum_Var_4                         = zeros(params.K, params.T, 'single');
  phi_E_4                           = cell(params.E, 1);
  phi_Var_4                         = cell(params.E, 1);

  for e = 1:params.E
    Ad                              = params.expByNumTasks(e);
    taskIndices                    = squeeze(params.expByTask(e,:));
    % N_{.tc.}
    % for E_q(N_{.tc.})
    phi_E_4{e}                      = bsxfun(@times, params.phi{e}, params.Fc{e});
    e_tmp_1                         = squeeze(sum(phi_E_4{e}, 1));
    if Ad == 1
        e_tmp_1                     = e_tmp_1';
    end;
    sum_E_4(:, taskIndices) = sum_E_4(:, taskIndices) + e_tmp_1;
    clear e_tmp_1;
    % for Var_q(N_{.tc.})
    phi_Var_4{e}                    = params.phi{e} .* (1 - params.phi{e});
    var_tmp_1                       = bsxfun(@times, phi_Var_4{e}, params.Fc{e});
    var_tmp_2                       = squeeze(sum(var_tmp_1, 1));
    if Ad == 1
      var_tmp_2 = var_tmp_2';
    end;
    sum_Var_4(:, taskIndices)    = sum_Var_4(:, taskIndices) + var_tmp_2;
    clear var_tmp_1 var_tmp_2
  end;

  if varBoundFlag
    % Compute the variational bound
    for t = 1:params.T
      n_4                             = single([0:(uint32(params.Nt(t))-1)]);
      log_tmp_4                       = log(single(params.alpha + [1:params.Nt(t)] - 1));
      for c = 1:params.K
        mu_4                          = single(sum_E_4(c, t));
        sigma_4                       = single(sqrt(sum_Var_4(c, t)));

        if sigma_4 <= 0
          newVarBound                 = newVarBound + sum(log_tmp_4(1:floor(mu_4)));
        else
          q_4                           = 1 - normcdf(n_4, mu_4, sigma_4);
          newVarBound                 = newVarBound + sum(q_4 .* log_tmp_4);
        end;

        clear mu_4 sigma_4;
      end
      clear n_4 log_tmp_4
    end;
  end;

  % Deduct variables of e, n
  for e = 1:params.E
    Nd                              = params.fociCounts(e);
    taskIndices                    = squeeze(params.expByTask(e,:));
    % N_{.tc.}
    % for E_q(N_{.tc.})
    tmp_1                           = repmat(sum_E_4(:, taskIndices), [1 1 Nd]);
    tmp_2                           = permute(tmp_1, [3 1 2]);
    E_4{e} = tmp_2 - phi_E_4{e};
    clear tmp_1 tmp_2;
    % for Var_q(N_{.tc.})
    tmp_3                           = repmat(sum_Var_4(:, taskIndices), [1 1 Nd]);
    tmp_4                           = permute(tmp_3, [3 1 2]);
    Var_4{e}                        = tmp_4 - phi_Var_4{e};
    % clear
    clear tmp_3 tmp_4
  end;
  clear sum_E_4 sum_Var_4 phi_E_4 phi_Var_4

  % update phi
  for e = 1:params.E
    X_4_tmp                         = params.alpha + E_4{e};
    X_4                             = log(X_4_tmp);
    Y_4                             = Var_4{e} ./ (X_4_tmp.^2);

    newPhi{e}                      = newPhi{e} + X_4;
    newPhi{e}                      = newPhi{e} - 0.5 * Y_4;

    clear X_4 X_4_tmp Y_4
  end;
  clear E_4 Var_4
