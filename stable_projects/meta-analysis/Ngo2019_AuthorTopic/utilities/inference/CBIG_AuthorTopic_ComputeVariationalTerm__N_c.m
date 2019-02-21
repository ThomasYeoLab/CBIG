function [newPhi, newVarBound] = ...
    CBIG_AuthorTopic_ComputeVariationalTerm__N_c(params, newPhi, ...
    newVarBound, varBoundFlag)
% [newPhi, newVarBound] = CBIG_AuthorTopic_ComputeVariationalTerm__N_c(
%    params, newPhi, newVarBound, varBoundFlag)
%
% Compute the expectation and variance of the count N_{..c.}^{-ef} in
% the update equation for the parameters of the variation distribution phi
% Optionally update the lower bound of the data log likelihood.
%
% Input:
%  - params      : struct containing parameters of the author-topic model
%                  and the Collapsed Variational Bayes (CVB) algorithm.
%  - newPhi      : current estimate of the parameters of the variational
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
%     CBIG_AuthorTopic_ComputeVariationalTerm__N_c(params, newPhi, ...
%     newVarBound, true)
% Update the parameters of the variational distribution and the
% variational bound of the data log likelihood by computing the
% expectation and variance of the count N_{..c.}^{-ef}.
%
% Reference:
%   1) Eq(14) of the Supplemental S1 of Ngo et al., 2019,
%   "Beyond Consensus: Embracing Heterogeneity in Curated
%   Neuroimaging Meta-Analysis"
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_1               = cell(params.E, 1);
  Var_1             = cell(params.E, 1);
  sum_E_1           = zeros(1, params.K, 'single');
  sum_Var_1         = zeros(1, params.K, 'single');
  phi_E_1           = cell(params.E, 1);
  phi_Var_1         = cell(params.E, 1);


  for e = 1:params.E
    %N_{..k.}
    % for E_q(N_{..k.})
    phi_E_1{e}      = squeeze(sum(params.phi{e},3));
    e_tmp_1         = bsxfun(@times, phi_E_1{e}, params.Fc{e});
    e_tmp_2         = squeeze(sum(e_tmp_1, 1));
    sum_E_1         = sum_E_1 + e_tmp_2;
    clear e_tmp_1 e_tmp2;
    % for Var_q(N_{..k.})
    phi_Var_1{e}    = phi_E_1{e} .* (1 - phi_E_1{e});
    var_tmp_1       = bsxfun(@times, phi_Var_1{e}, params.Fc{e});
    var_tmp_2       = squeeze(sum(var_tmp_1, 1));
    sum_Var_1       = sum_Var_1 + var_tmp_2;
    % clear
    clear var_tmp_1 var_tmp_2;
  end;

  % Compute the variational bound
  if varBoundFlag
    n_1                 = single([0:(params.totalFociCount-1)]);
    log_tmp_1           = log(single(params.V * params.eta + [1:params.totalFociCount] - 1));
    for k = 1:params.K
      mu_1              = single(sum_E_1(k));
      sigma_1           = single(sqrt(squeeze(sum_Var_1(k))));

      if sigma_1 <= 0
        newVarBound   = newVarBound - sum(log_tmp_1(1:floor(mu_1)));
      else
        q_1               = 1 - normcdf(n_1, mu_1, sigma_1);
        newVarBound     = newVarBound - sum(q_1 .* log_tmp_1);
      end;

      clear mu_1 sigma_1 q_1;
    end;
    clear n_1 log_tmp_1;
  end;

  %.Eeduct variables of d, n
  for e = 1:params.E
    Nd              = params.fociCounts(e);
    %N_{..k.}
    % for E_q(N_{..k.})
    tmp_1           = repmat(sum_E_1, [Nd 1]);
    E_1{e}          = tmp_1 - phi_E_1{e};
    clear tmp_1;
    % for Var_q(N_{..k.})
    tmp_2           = repmat(sum_Var_1, [Nd 1]);
    Var_1{e}        = tmp_2 - phi_Var_1{e};
    % clear
    clear tmp_2;
  end;
  clear sum_E_1 sum_Var_1 phi_E_1 phi_Var_1

  % update phi
  for e = 1:params.E
    X_1_tmp         = params.V * params.eta + E_1{e};
    X_1             = log(X_1_tmp);

    Y_1             = Var_1{e} ./ (X_1_tmp.^2);

    newPhi{e}      = bsxfun(@minus, newPhi{e}, X_1);
    newPhi{e}      = bsxfun(@plus, newPhi{e}, 0.5 * Y_1);

    clear X_1_tmp X_1 Y_1
  end;
  %clear
  clear E_1 Var_1
