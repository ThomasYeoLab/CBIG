function [newPhi, newVarBound] = ...
  CBIG_AuthorTopic_ComputeVariationalTerm__N_cv(params, newPhi, ...
  newVarBound, varBoundFlag)
% [newPhi, newVarBound] = ...
%   CBIG_AuthorTopic_ComputeVariationalTerm__N_cv(params, newPhi, ...
%   newVarBound, varBoundFlag)
%
% Compute the expectation and variance of the the count
% N_{..cv_{ef}}^{-ef} in the update equation for the parameter of the
% variation distribution phi.
% Optionally update the variational bound of the data log likelihood.
%
% Input:
%  - params      : struct containing parameters of the author-topic model
%                  and the Collapsed Variational Bayes (CVB) algorithm.
%  - newPhi      : current estimate of the paramters of the variational
%                  distribution.
%  - newVarBound : current estimate of the lower bound of the data log
%                  likelihood.
%  - varBoundFlag: if true, update the lower bound and do nothing
%                  otherwise.
% Output:
%  - newPhi: updated estimate of the parameters of the variational
%            distribution.
%  - newVarBound: updated estimate of the variational bound.
%
% Example:
%   [newPhi, newVarBound] = ...
%     CBIG_AuthorTopic_ComputeVariationalTerm__N_cv(params, newPhi, ...
%     newVarBound, true)
% Update the parameters of the variational distribution and the
% variational bound of the data log likelihood by computing the
% expectation and variance of the count N_{..cv_{ef}}^{-ef}.
%
% Reference:
%   1) Eq(14) of the Supplemental S1 of Ngo et al., 2019,
%   "Beyond Consensus: Embracing Heterogeneity in Neuroimaging
%    Meta-Analysis"
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


  E_2                           = cell(params.E, 1);
  Var_2                         = cell(params.E, 1);
  sum_E_2                       = zeros(params.V, params.K, 'single');
  sum_Var_2                     = zeros(params.V, params.K, 'single');
  phi_E_2                       = cell(params.E, 1);
  phi_Var_2                     = cell(params.E, 1);

  for e = 1:params.E
    fociIndices                 = params.fociIndices{e};
    % N_{..cv_{ef}}
    % for E_q(N_{..cv_{ef}})
    phi_E_2{e}                  = squeeze(sum(params.phi{e},3));
    e_tmp_1                     = bsxfun(@times, phi_E_2{e}, params.Fc{e});
    sum_E_2(fociIndices,:)      = sum_E_2(fociIndices,:) + e_tmp_1;
    clear e_tmp_1;
    % for Var_q(N_{..cv_{ef}})
    phi_Var_2{e}                = phi_E_2{e} .* (1 - phi_E_2{e});
    var_tmp_1                   = bsxfun(@times, phi_Var_2{e}, params.Fc{e});
    sum_Var_2(fociIndices,:)    = sum_Var_2(fociIndices,:) + var_tmp_1;
    % clear
    clear var_tmp_1;
  end;

  if varBoundFlag
    % Compute the variational bound
    for i = 1:size(params.uniqueNv, 1)
      word_count                  = params.uniqueNv(i);
      if word_count >= 1
        fociIndices               = params.fociIndicesByCount{i};

        e_tmp                     = single(sum_E_2(fociIndices,:));
        e_tmp                     = e_tmp(:);
        var_tmp                   = single(sum_Var_2(fociIndices,:));
        var_tmp                   = var_tmp(:);

        if size(e_tmp, 1) == 1
            e_tmp = e_tmp';
        end;
        if size(var_tmp, 1) == 1
            var_tmp = var_tmp';
        end;

        positiveVarIndices        = find(var_tmp > 0);
        zeroVarIndices            = find(var_tmp == 0);
        mu_2                      = repmat(e_tmp, [1 word_count]);
        sigma_2                   = repmat(sqrt(var_tmp), [1 word_count]);
        clear e_tmp var_tmp;

        n_2                       = repmat(single(0:(word_count-1)), [size(mu_2, 1) 1]);
        log_tmp_2                 = log(params.eta + single(1:word_count) - 1);

        q_2                       = zeros(size(mu_2), 'single');
        q_2(positiveVarIndices, :)   = 1 - normcdf(n_2(positiveVarIndices,:), mu_2(positiveVarIndices,:), sigma_2(positiveVarIndices,:));
        tmp                       = bsxfun(@times, q_2, log_tmp_2);
        newVarBound               = newVarBound + sum(tmp(:));
        clear tmp;

        for j = 1:numel(zeroVarIndices)
          idx                     = zeroVarIndices(j);
          newVarBound             = newVarBound + sum(log_tmp_2(1:floor(mu_2(idx,1))));
        end;
      end;
    end;
    clear mu_2 sigma_2 q_2 n_2 log_tmp_2;
  end;

  %.Eeduct variables of d, n
  for e = 1:params.E
    fociIndices                 = params.fociIndices{e};
    % N_{..cv_{ef}}
    % for E_q(N_{..cv_{ef}})
    E_2{e}                      = sum_E_2(fociIndices, :) - phi_E_2{e};
    % for Var_q(N_{..cv{ef}})
    Var_2{e}                    = sum_Var_2(fociIndices, :) - phi_Var_2{e};
  end;
  clear sum_E_2 sum_Var_2 phi_E_2 phi_Var_2;

  % update phi
  for e = 1:params.E
    X_2_tmp                     = params.eta + E_2{e};
    X_2                         = log(X_2_tmp);
    Y_2                         = Var_2{e} ./ (X_2_tmp.^2);

    newPhi{e}                   = bsxfun(@plus, newPhi{e}, X_2);
    newPhi{e}                   = bsxfun(@minus, newPhi{e}, 0.5 * Y_2);

    clear X_2_tmp X_2 Y_2
  end;
  %clear
  clear E_2 Var_2
