function params = CreateEmptyAT_CVBParams(SEED, K, BASE_DIR, DATA_PATH, BURN_IN)
  % params = CreateEmptyAT_CVBParams(SEED, K, BASE_DIR, DATA_PATH, BURN_IN)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  params.mask_path          = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'meta-analysis', 'Ngo2017_AuthorTopic', 'utilities', 'brain_mask', 'MNI_mask_conformed.2mm.0.1.nii.gz');
  params.base_dir           = BASE_DIR;
  params.data_path          = DATA_PATH;
  params.em_min_iter        = 5;
  params.e_step_max_iter    = 1;
  
  if nargin < 5
    params.var_bound_burn_in = 26;
    params.em_max_iter       = 100;
  else
    params.var_bound_burn_in = BURN_IN;
    params.em_max_iter       = BURN_IN * 2;
    if params.em_max_iter < 100
      params.em_max_iter  = 100;
    end
  end

  params.var_bound_interval = 5;
  params.var_bound_convergence  = 5e-4;

  params.lower_bound    = -inf;

  params.seed = SEED;
  params.K = K;
  
  params.smooth_kernel_size = 31;
  params.smooth_kernel_sigma = 11.0;
  params.smooth_repeats = 1;
