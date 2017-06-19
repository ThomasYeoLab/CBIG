function [new_phi, new_var_bound] = ComputeThirdTerm__N_a(params, new_phi, new_var_bound, var_bound_flag)
  % [new_phi, new_var_bound] = ComputeThirdTerm__N_a(params, new_phi, new_var_bound, var_bound_flag)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_3                           = cell(params.D, 1);
  Var_3                         = cell(params.D, 1);
  sum_E_3                       = zeros(1, params.A, 'single');
  sum_Var_3                     = zeros(1, params.A, 'single');
  phi_E_3                       = cell(params.D, 1);
  phi_Var_3                     = cell(params.D, 1);

  for d = 1:params.D
    author_indices              = squeeze(params.author_indices(d,:));
    % N_{.a..}
    % for E_q(N_{.a..})
    phi_E_3{d}                  = squeeze(sum(params.phi{d},2));
    e_tmp_1                     = bsxfun(@times, phi_E_3{d}, params.Wc{d});
    e_tmp_2                     = squeeze(sum(e_tmp_1, 1));
    sum_E_3(author_indices)     = sum_E_3(author_indices) + e_tmp_2;
    clear e_tmp_1 e_tmp_2;
    % for Var_q(N_{kw_{dn}})
    phi_Var_3{d}                = phi_E_3{d} .* (1 - phi_E_3{d});
    var_tmp_1                   = bsxfun(@times, phi_Var_3{d}, params.Wc{d});
    var_tmp_2                   = squeeze(sum(var_tmp_1, 1));
    sum_Var_3(author_indices)   = sum_Var_3(author_indices) + var_tmp_2;
    % clear
    clear var_tmp_1 var_tmp_2;
  end;

  if var_bound_flag
    % Compute the variational bound
    for a = 1:params.A
      n_3               = single([0:(uint32(params.Na(a)) - 1)]);  
      log_tmp_3         = log(single(params.K * params.alpha + [1:params.Na(a)] - 1));
      mu_3              = single(sum_E_3(a));
      sigma_3           = single(sqrt(sum_Var_3(a)));
      
      if sigma_3 <= 0
        new_var_bound   = new_var_bound - sum(log_tmp_3(1:floor(mu_3)));
      else
        q_3             = 1 - normcdf(n_3, mu_3, sigma_3);
        new_var_bound   = new_var_bound - sum(q_3 .* log_tmp_3);
      end;

      clear mu_3 sigma_3;
    end;
    clear n_3 log_tmp_3
  end;
  

  % Deduct variables of d, n
  for d = 1:params.D
    Nd                          = params.word_counts(d);
    author_indices              = squeeze(params.author_indices(d,:));
    % N_{.a..}
    % for E_q(N_{.a..})
    tmp_1                       = repmat(sum_E_3(author_indices), [Nd 1]);
    E_3{d}                      = tmp_1 - phi_E_3{d};
    clear tmp_1;
    % for Var_q(N_{.a..})
    tmp_2                       = repmat(sum_Var_3(author_indices), [Nd 1]);
    Var_3{d}                    = tmp_2 - phi_Var_3{d};
    % clear
    clear tmp_2
  end;
  clear sum_E_3 sum_Var_3 phi_E_3 phi_Var_3

  % update phi
  for d = 1:params.D  
    X_3_tmp                     = params.K * params.alpha + E_3{d};
    X_3                         = log(X_3_tmp);
    X_3                         = permute(X_3, [1 3 2]);
    Y_3_tmp                     = Var_3{d} ./ (X_3_tmp.^2);
    Y_3                         = permute(Y_3_tmp, [1 3 2]);

    new_phi{d}                  = bsxfun(@minus, new_phi{d}, X_3);
    new_phi{d}                  = bsxfun(@plus, new_phi{d}, 0.5 * Y_3);

    clear X_3 X_3_tmp Y_3 Y_3_tmp
  end;
  clear E_3 Var_3
