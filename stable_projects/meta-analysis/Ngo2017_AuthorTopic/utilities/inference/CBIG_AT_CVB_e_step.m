function [new_phi, new_var_bound] = CBIG_AT_CVB_e_step(params, new_var_bound, var_bound_flag)
  % [new_phi, new_var_bound] = CBIG_AT_CVB_e_step(params, new_var_bound, var_bound_flag)
  % params.Wb = sparse D x V logical array, with Wb(d,v)  = 1 
  % if the v-th vocabulary word is in document d
  %
  % params.Wc = cell array of size D x 1, where Wc{d} = Nd x 1 vector
  % with Wc{d}(n) = number of times the n-th unique word appears in the
  % d-th document
  %
  % params.doc_by_author = sparse D x A logical array, with
  % doc_by_author(d,a) = 1 if author a is in document d

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  new_phi = cell(params.D, 1);  
  for iter = 1:params.e_step_max_iter
    for d = 1:params.D
      Nd = params.word_counts(d);
      Ad = params.doc_by_num_authors(d);
      new_phi{d} = zeros(Nd, params.K, Ad, 'single');
    end;

    [new_phi, new_var_bound] = CBIG_ComputeFirstTerm__N_k(params, new_phi, new_var_bound, var_bound_flag);
    [new_phi, new_var_bound] = CBIG_ComputeSecondTerm__N_kw(params, new_phi, new_var_bound, var_bound_flag);
    [new_phi, new_var_bound] = CBIG_ComputeThirdTerm__N_a(params, new_phi, new_var_bound, var_bound_flag);
    [new_phi, new_var_bound] = CBIG_ComputeFourthTerm__N_ak(params, new_phi, new_var_bound, var_bound_flag);

    % normalize new phi
    new_phi = CBIG_normalizePhi(new_phi, params);
  end;
