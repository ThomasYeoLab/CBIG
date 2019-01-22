function CBIG_AuthorTopic_RunInference(seed, K, alpha, eta, doSmoothing, baseDirectory, dataPath, burnin, isSavingPhi)
  % CBIG_AuthorTopic_RunInference(seed, K, alpha, eta, doSmoothing, baseDirectory, dataPath, burnin, isSavingPhi)
  %
  % Wrapper function to estimate the author-topic models parameters using the Collapsed Variational Bayes (CVB) algorithm.
  %
  % Input:
  %  - seed: random initialization (seed) used for the current run of the CVB algorithm.
  %  - K: number of components of the author-topic model.
  %  - alpha & eta: hyperparameters of the Dirichlet priors of the author-topic model
  %  - doSmoothing: if doSmoothing > 1, performs spatial smoothing of the CVB algorithm's variational distribution. If doSmoothing = 0, do not perform spatial smoothing.
  %  - baseDirectory: absolute path to the base directory containing the output of the current CVB algorithm's run.
  %  - dataPath: absolute path to the .mat file containing the author-topic model's input data.
  %  - burnin: number of iterations performed by the CVB algorithm before checking for convergence. The default is set by CBIG_AuthorTopic_SetupParameters.m
  %  - isSavingPhi: if isSavingPhi > 1, the variational distribution of the CVB algorithm is saved. If isSavingPhi = 0, the variational distribution of the CVB algorithm is not saved (e.g. to save storage space). Default value: 0.
  %  Output:
  %  - The estimates of the author-topic model's parameters are saved at <baseDirectory>/outputs/K<K>/alpha<alpha>_eta<eta>/params_K<K>_SEED<seed>.mat
  %
  % Example:
  %   CBIG_AuthorTopic_RunInference(23, 2, 100, 0.01, 1, '/AT_outputs', '/data/selfGeneratedData.mat', 50, 1)
  %   Estimate the author-topic model's parameters by the CVB algorithm. The CVB algorithm is initialized with the random seed 23, hyperparameters alpha = 100 & eta = 0.01. The author-topic model is assumed to have 2 components. The CVB algorithm's variational distribution is spatially smoothed and has its output saved at /AT_outputs. The algorithm uses input saved in /data/selfGeneratedData.mat. The CVB algorithm runs for 50 iterations before checking for convergence and its variational distribution is saved in the output. The output file is /AT_outputs/outputs/K2/alpha100_eta0.01/params_K2_SEED23.mat
  %
  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % Format arguments if necessary
  if (ischar(seed))
    seed = str2num(seed); %#ok<*ST2NM>
  end;

  if (ischar(K))
    K = str2num(K);
  end;

  if (ischar(alpha))
    alpha = str2num(alpha);
  end;

  if (ischar(eta))
    eta = str2num(eta);
  end;

  if (ischar(doSmoothing))
    doSmoothing = str2num(doSmoothing);
  end

  if nargin == 8 && (ischar(burnin))
    burnin = str2num(burnin);
  end

  if nargin == 9 && (ischar(isSavingPhi))
    isSavingPhi = str2num(isSavingPhi);
  else
    isSavingPhi = 0;
  end

  % set up the model parameters
  if nargin < 8
    params                = CBIG_AuthorTopic_SetupParameters(seed, K, baseDirectory, dataPath);
  else
    params                = CBIG_AuthorTopic_SetupParameters(seed, K, baseDirectory, dataPath, burnin);
  end

  rng(params.seed * 10000, 'twister');
  params                = CBIG_AuthorTopic_InitializeParams(params, alpha, eta);

  % initial smoothing
  if doSmoothing ~= 0
    disp('Smoothing variational distribution');

    params.phi = CBIG_AuthorTopic_SmoothPhi(params);
  end


  allVarBounds        = zeros(params.maxIteration, 1, 'single');

  for iter = 1:params.maxIteration
    disp(['Iteration ' num2str(iter)]);

    % initialize new variational bound
    newVarBound             = 0;

    % for efficiency, only compute variational bound every varBoundInterval iterations
    varBoundFlag            = false;
    if (rem(iter, params.varBoundInterval) == 0)
      varBoundFlag            = true;
      for e = 1:params.E
        tmp                     = params.phi{e} .* log(params.phi{e});
        tmp(params.phi{e} == 0) = 0;
        newVarBound           = newVarBound - sum(sum(sum(tmp, 3), 2));
      end;
    end;

    [newPhi, newVarBound]  = CBIG_AuthorTopic_EstimateVariationalDistribution(params, newVarBound, varBoundFlag);
    allVarBounds(iter)      = newVarBound;

    % update the variational distribution
    params.phi                = newPhi;

    % check for convergence
    tmpVarbound = allVarBounds;
    tmpVarbound(tmpVarbound == 0) = [];
    if iter > params.varBoundBurnIn
      lastVarbound = tmpVarbound(numel(tmpVarbound));
      secondLastVarbound = tmpVarbound(numel(tmpVarbound) - 1);
      percentChange = abs((lastVarbound - secondLastVarbound) / secondLastVarbound);
      disp(['Variational bound: ' num2str(lastVarbound)]);
      disp(['Percentage change of variational bound: ' num2str(percentChange)]);
      if percentChange < params.varBoundConvergenceThresh
        disp(['Convergence occured in iteration: ' num2str(iter)]);
        break;
      end;
    end;
    clear tmpVarbound lastVarbound secondLastVarbound percentChange

    params.iter = iter;
  end;

  % estimate model parameters
  [theta, beta]         = CBIG_AuthorTopic_EstimateParams(params);
  params.theta          = theta;
  params.beta           = beta;
  clear theta beta;

  allVarBounds(allVarBounds == 0) = [];
  params.varBound            = allVarBounds;

  if ~isSavingPhi
    params                = rmfield(params, 'phi');
  end

  outputDir = fullfile(params.baseDir, 'outputs', ['K' num2str(params.K)]);
  mkdir(outputDir);
  targDir =['alpha' num2str(params.alpha) '_eta' num2str(params.eta)];
  mkdir(fullfile(outputDir, targDir));
  paramsFilePath            = fullfile(outputDir, targDir, ['params_K' num2str(params.K, '%3d') '_SEED' num2str(params.seed, '%3d') '.mat']);
  disp(['Saving model estimate to' paramsFilePath]);
  save(paramsFilePath, 'params');
