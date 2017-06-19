function [new_phi, new_var_bound] = ComputeFourthTerm__N_ak(params, new_phi, new_var_bound, var_bound_flag)
  % [new_phi, new_var_bound] = ComputeFourthTerm__N_ak(params, new_phi, new_var_bound, var_bound_flag)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_4                               = cell(params.D, 1);
  Var_4                             = cell(params.D, 1);
  sum_E_4                           = zeros(params.K, params.A, 'single');
  sum_Var_4                         = zeros(params.K, params.A, 'single');
  phi_E_4                           = cell(params.D, 1);
  phi_Var_4                         = cell(params.D, 1);

  for d = 1:params.D
    Ad                              = params.doc_by_num_authors(d);
    author_indices                  = squeeze(params.author_indices(d,:));
    % N_{.ak.}
    % for E_q(N_{.ak.})
    phi_E_4{d}                      = bsxfun(@times, params.phi{d}, params.Wc{d});        
    e_tmp_1                         = squeeze(sum(phi_E_4{d}, 1));
    if Ad == 1
        e_tmp_1                     = e_tmp_1';
    end;
    sum_E_4(:, author_indices) = sum_E_4(:, author_indices) + e_tmp_1;
    clear e_tmp_1;
    % for Var_q(N_{.ak.})
    phi_Var_4{d}                    = params.phi{d} .* (1 - params.phi{d});
    var_tmp_1                       = bsxfun(@times, phi_Var_4{d}, params.Wc{d});
    var_tmp_2                       = squeeze(sum(var_tmp_1, 1));
    if Ad == 1
        var_tmp_2 = var_tmp_2';
    end;
    sum_Var_4(:, author_indices)    = sum_Var_4(:, author_indices) + var_tmp_2;
    clear var_tmp_1 var_tmp_2
  end;

  if var_bound_flag
    % Compute the variational bound
    for a = 1:params.A
      n_4                             = single([0:(uint32(params.Na(a))-1)]);  
      log_tmp_4                       = log(single(params.alpha + [1:params.Na(a)] - 1));
      for k = 1:params.K
        mu_4                          = single(sum_E_4(k, a));
        sigma_4                       = single(sqrt(sum_Var_4(k, a)));
        
        if sigma_4 <= 0
          new_var_bound               = new_var_bound + sum(log_tmp_4(1:floor(mu_4)));
        else
          q_4                           = 1 - normcdf(n_4, mu_4, sigma_4);
          new_var_bound                 = new_var_bound + sum(q_4 .* log_tmp_4);
        end;

        clear mu_4 sigma_4;
      end
      clear n_4 log_tmp_4
    end;
  end;

  % Deduct variables of d, n
  for d = 1:params.D
    Nd                              = params.word_counts(d);
    author_indices                  = squeeze(params.author_indices(d,:));
    % N_{.ak.}
    % for E_q(N_{.ak.})
    tmp_1                           = repmat(sum_E_4(:,author_indices), [1 1 Nd]);
    tmp_2                           = permute(tmp_1, [3 1 2]);
    E_4{d} = tmp_2 - phi_E_4{d};
    clear tmp_1 tmp_2;
    % for Var_q(N_{.ak.})
    tmp_3                           = repmat(sum_Var_4(:,author_indices), [1 1 Nd]);
    tmp_4                           = permute(tmp_3, [3 1 2]);
    Var_4{d}                        = tmp_4 - phi_Var_4{d};
    % clear
    clear tmp_3 tmp_4
  end;
  clear sum_E_4 sum_Var_4 phi_E_4 phi_Var_4

  % update phi
  for d = 1:params.D
    X_4_tmp                         = params.alpha + E_4{d};
    X_4                             = log(X_4_tmp);
    Y_4                             = Var_4{d} ./ (X_4_tmp.^2);

    new_phi{d}                      = new_phi{d} + X_4;
    new_phi{d}                      = new_phi{d} - 0.5 * Y_4;

    clear X_4 X_4_tmp Y_4
  end;
  clear E_4 Var_4


