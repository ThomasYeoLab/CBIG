function RunAT_CVB(SEED, K, ALPHA, ETA, SMOOTH, BASE_DIR, DATA_PATH, BURN_IN, SAVE_PHI)
  % RunAT_CVB(SEED, K, ALPHA, ETA, SMOOTH, BASE_DIR, DATA_PATH, BURN_IN, SAVE_PHI)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
  tic;

  if (ischar(SEED))
    SEED = str2num(SEED);
  end;
  
  if (ischar(K))
    K = str2num(K);
  end;

  if (ischar(ALPHA))
    ALPHA = str2num(ALPHA);
  end;

  if (ischar(ETA))
    ETA = str2num(ETA);
  end;

  if (ischar(SMOOTH))
    SMOOTH = str2num(SMOOTH);
  end
  
  if nargin == 8 && (ischar(BURN_IN))
    BURN_IN = str2num(BURN_IN);
  end
  
  if nargin == 9 && (ischar(SAVE_PHI))
    SAVE_PHI = str2num(SAVE_PHI);
  else
    SAVE_PHI = 0;
  end

  % initialize params
  if nargin < 8
    params                = CBIG_CreateEmptyAT_CVBParams(SEED, K, BASE_DIR, DATA_PATH);
  else
    params                = CBIG_CreateEmptyAT_CVBParams(SEED, K, BASE_DIR, DATA_PATH, BURN_IN);
  end

  rng(sum(10000*params.seed), 'twister');
  params                = CBIG_InitAT_CVBParams(params, ALPHA, ETA);
  
  % initial smoothing
  if SMOOTH ~= 0
    disp('Volume smoothing');
    % params.phi = SmoothPhi(params);
    params.phi = CBIG_SmoothPhiWithPredefinedFFTGaussianKernel(params);
    % params.phi = SmoothPhiSameNumAuthorPerDoc(params);
  end

  
  all_var_bounds        = zeros(params.em_max_iter, 1, 'single');

  for iter = 1:params.em_max_iter
    % e-step
    disp('');
    disp(['EM Iteration ' num2str(iter)]);

    % initialize new variational bound
    new_var_bound             = 0; 

    % only compute variational bound every VAR_BOUND_INTERVAL iterations
    var_bound_flag            = false;
    if (rem(iter, params.var_bound_interval) == 0)
      var_bound_flag            = true;
      for d = 1:params.D
        tmp                     = params.phi{d} .* log(params.phi{d});
        tmp(params.phi{d} == 0) = 0;
        new_var_bound           = new_var_bound - sum(sum(sum(tmp, 3), 2));
      end;
    end;
    
    [new_phi, new_var_bound]  = CBIG_AT_CVB_e_step(params, new_var_bound, var_bound_flag);
    all_var_bounds(iter)      = new_var_bound;

    % update phi
    params.phi                = new_phi;

    % Check for convergence
    tmpVarbound = all_var_bounds;
    tmpVarbound(tmpVarbound == 0) = [];
    if iter > params.var_bound_burn_in
      lastVarbound = tmpVarbound(numel(tmpVarbound));
      secondLastVarbound = tmpVarbound(numel(tmpVarbound) - 1);
      percentChange = abs((lastVarbound - secondLastVarbound) / secondLastVarbound);
      disp(['Variational bound: ' num2str(lastVarbound)]);
      disp(['Percentage change of variational bound: ' num2str(percentChange)]);
      if percentChange < params.var_bound_convergence
        disp(['Convergence occured in iteration: ' num2str(iter)]);
        break;
      end;
    end;
    clear tmpVarbound lastVarbound secondLastVarbound percentChange

    params.iter = iter;
  end;

  % estimate initial parameters
  [theta, beta]         = CBIG_EstimateATParams(params);
  params.theta          = theta;
  params.beta           = beta;
  clear theta beta;

  all_var_bounds(all_var_bounds == 0) = [];
  params.var_bound            = all_var_bounds;
  params.elapsedTime          = toc;
  
  if ~SAVE_PHI
    params                = rmfield(params, 'phi');
  end
 
  outputDir = fullfile(params.base_dir, 'outputs', ['K' num2str(params.K)]);
  mkdir(outputDir);
  targDir =['alpha' num2str(params.alpha) '_eta' num2str(params.eta)];
  mkdir(fullfile(outputDir, targDir));
  params_file_path            = fullfile(outputDir, targDir, ['params_K' num2str(params.K, '%3d') '_SEED' num2str(params.seed, '%3d') '.mat']);
  disp(params_file_path);    
  save(params_file_path, 'params');
