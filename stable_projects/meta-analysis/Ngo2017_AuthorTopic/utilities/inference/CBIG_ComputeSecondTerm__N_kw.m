function [new_phi, new_var_bound] = ComputeSecondTerm__N_kw(params, new_phi, new_var_bound, var_bound_flag)
  % [new_phi, new_var_bound] = ComputeSecondTerm__N_kw(params, new_phi, new_var_bound, var_bound_flag)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  E_2                           = cell(params.D, 1);
  Var_2                         = cell(params.D, 1);
  sum_E_2                       = zeros(params.V, params.K, 'single');
  sum_Var_2                     = zeros(params.V, params.K, 'single');
  phi_E_2                       = cell(params.D, 1);
  phi_Var_2                     = cell(params.D, 1);

  for d = 1:params.D
    word_indices                = params.word_indices{d};
    % N_{..kw_{dn}}
    % for E_q(N_{..kw_{dn}})
    phi_E_2{d}                  = squeeze(sum(params.phi{d},3));
    e_tmp_1                     = bsxfun(@times, phi_E_2{d}, params.Wc{d});
    sum_E_2(word_indices,:)     = sum_E_2(word_indices,:) + e_tmp_1;
    clear e_tmp_1;
    % for Var_q(N_{kw_{dn}})
    phi_Var_2{d}                = phi_E_2{d} .* (1 - phi_E_2{d});
    var_tmp_1                   = bsxfun(@times, phi_Var_2{d}, params.Wc{d});
    sum_Var_2(word_indices,:)   = sum_Var_2(word_indices,:) + var_tmp_1;
    % clear
    clear var_tmp_1;
  end;
  
  if var_bound_flag
    % Compute the variational bound
    for i = 1:size(params.unique_Nv, 1)
      word_count                  = params.unique_Nv(i);
      if word_count >= 1
        word_indices              = params.word_indices_by_count{i};
        
        e_tmp                     = single(sum_E_2(word_indices,:));
        e_tmp                     = e_tmp(:);
        var_tmp                   = single(sum_Var_2(word_indices,:));
        var_tmp                   = var_tmp(:);
        
        if size(e_tmp, 1) == 1
            e_tmp = e_tmp';
        end;
        if size(var_tmp, 1) == 1
            var_tmp = var_tmp';
        end;

        pos_var_indices           = find(var_tmp > 0);
        zero_var_indices          = find(var_tmp == 0);
        mu_2                      = repmat(e_tmp, [1 word_count]);
        sigma_2                   = repmat(sqrt(var_tmp), [1 word_count]);
        clear e_tmp var_tmp;
        
        n_2                       = repmat(single(0:(word_count-1)), [size(mu_2, 1) 1]);  
        log_tmp_2                 = log(params.eta + single(1:word_count) - 1);

        q_2                       = zeros(size(mu_2), 'single');
        q_2(pos_var_indices, :)   = 1 - normcdf(n_2(pos_var_indices,:), mu_2(pos_var_indices,:), sigma_2(pos_var_indices,:));
        tmp                       = bsxfun(@times, q_2, log_tmp_2);
        new_var_bound             = new_var_bound + sum(tmp(:));
        clear tmp;
        
        for j = 1:numel(zero_var_indices)
          idx                     = zero_var_indices(j);
          new_var_bound           = new_var_bound + sum(log_tmp_2(1:floor(mu_2(idx,1))));
        end;
      end;
    end;
    clear mu_2 sigma_2 q_2 n_2 log_tmp_2;
  end;
  
  % Deduct variables of d, n
  for d = 1:params.D
    word_indices                = params.word_indices{d};
    % N_{..kw_{dn}}
    % for E_q(N_{..kw_{dn}})
    E_2{d}                      = sum_E_2(word_indices, :) - phi_E_2{d};
    % for Var_q(N_{..kw{dn}})
    Var_2{d}                    = sum_Var_2(word_indices, :) - phi_Var_2{d};
  end;
  clear sum_E_2 sum_Var_2 phi_E_2 phi_Var_2;

  % update phi
  for d = 1:params.D
    X_2_tmp                     = params.eta + E_2{d};
    X_2                         = log(X_2_tmp);
    Y_2                         = Var_2{d} ./ (X_2_tmp.^2);

    new_phi{d}                  = bsxfun(@plus, new_phi{d}, X_2);
    new_phi{d}                  = bsxfun(@minus, new_phi{d}, 0.5 * Y_2);

    clear X_2_tmp X_2 Y_2
  end;
  %clear
  clear E_2 Var_2
