function [new_phi, new_var_bound] = CBIG_ComputeFirstTerm__N_k(params, new_phi, new_var_bound, var_bound_flag)
  % [new_phi, new_var_bound] = CBIG_ComputeFirstTerm__N_k(params, new_phi, new_var_bound, var_bound_flag)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_1               = cell(params.D, 1);
  Var_1             = cell(params.D, 1);
  sum_E_1           = zeros(1, params.K, 'single');
  sum_Var_1         = zeros(1, params.K, 'single');
  phi_E_1           = cell(params.D, 1);
  phi_Var_1         = cell(params.D, 1);  
    

  for d = 1:params.D
    %N_{..k.}
    % for E_q(N_{..k.})
    phi_E_1{d}      = squeeze(sum(params.phi{d},3));
    e_tmp_1         = bsxfun(@times, phi_E_1{d}, params.Wc{d});
    e_tmp_2         = squeeze(sum(e_tmp_1, 1));
    sum_E_1         = sum_E_1 + e_tmp_2;
    clear e_tmp_1 e_tmp2;
    % for Var_q(N_{..k.})
    phi_Var_1{d}    = phi_E_1{d} .* (1 - phi_E_1{d});
    var_tmp_1       = bsxfun(@times, phi_Var_1{d}, params.Wc{d});
    var_tmp_2       = squeeze(sum(var_tmp_1, 1));
    sum_Var_1       = sum_Var_1 + var_tmp_2;
    % clear
    clear var_tmp_1 var_tmp_2;
  end;

  % Compute the variational bound
  if var_bound_flag
    n_1                 = single([0:(params.total_words_count-1)]);  
    log_tmp_1           = log(single(params.V * params.eta + [1:params.total_words_count] - 1));
    for k = 1:params.K
      mu_1              = single(sum_E_1(k));
      sigma_1           = single(sqrt(squeeze(sum_Var_1(k))));

      if sigma_1 <= 0
        new_var_bound   = new_var_bound - sum(log_tmp_1(1:floor(mu_1)));
      else
        q_1               = 1 - normcdf(n_1, mu_1, sigma_1);
        new_var_bound     = new_var_bound - sum(q_1 .* log_tmp_1);
      end;

      clear mu_1 sigma_1 q_1;
    end;
    clear n_1 log_tmp_1;
  end;
  
  % Deduct variables of d, n
  for d = 1:params.D
    Nd              = params.word_counts(d);    
    %N_{..k.}
    % for E_q(N_{..k.})
    tmp_1           = repmat(sum_E_1, [Nd 1]);
    E_1{d}          = tmp_1 - phi_E_1{d};
    clear tmp_1;
    % for Var_q(N_{..k.})
    tmp_2           = repmat(sum_Var_1, [Nd 1]);
    Var_1{d}        = tmp_2 - phi_Var_1{d};
    % clear
    clear tmp_2;
  end;
  clear sum_E_1 sum_Var_1 phi_E_1 phi_Var_1

  % update phi
  for d = 1:params.D
    X_1_tmp         = params.V * params.eta + E_1{d};
    X_1             = log(X_1_tmp);

    Y_1             = Var_1{d} ./ (X_1_tmp.^2);

    new_phi{d}      = bsxfun(@minus, new_phi{d}, X_1);
    new_phi{d}      = bsxfun(@plus, new_phi{d}, 0.5 * Y_1);

    clear X_1_tmp X_1 Y_1
  end;
  %clear         
  clear E_1 Var_1
