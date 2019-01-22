function params = CBIG_AuthorTopic_SetupParameters( ...
    seed, K, baseDirectory, dataPath, burnin)
% params = CBIG_AuthorTopic_SetupParameters( ...
%   seed, K, baseDirectory, dataPath, burnin)
%
% Set up the basic parameters of the Collapsed Variational Bayes (CVB)
% algorithm for the author-topic model.
%
% Input:
%  - seed         : seed of the random number generator used in the current
%                   run of the CVB algorithm.
%  - K            : number of components of the author-topic model.
%  - baseDirectory: absolute path to the base directory containing the
%                   output of the current CVB algorithm's run.
%  - dataPath     : absolute path to the .mat file containing the input
%                   meta-analytic data.
%  - burnin       : number of iterations performed by the CVB algorithm
%                   before checking for convergence. Default value: 26.
%
% Example:
%   params = CBIG_AuthorTopic_SetupParameters(23, 2, '/AT_outputs', ...
%              '/data/setGeneratedData.mat', 50)
%   Set up the CVB algorithm's parameters. The CVB algorithm is initialized
%   with the random seed 23. The author-topic model is assumed to have 2
%   components. The CVB algorithm has its output saved at /AT_outputs.
%   The algorithm uses input saved in /data/selfGeneratedData.mat. The CVB
%   algorithm runs for 50 iterations before checking for convergence.
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  params.maskPath          = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', ...
  'meta-analysis', 'Ngo2019_AuthorTopic', 'utilities', 'mask', 'MNI_mask_conformed.2mm.0.1.nii.gz');
  params.baseDir           = baseDirectory;
  params.dataPath          = dataPath;

  params.seed = seed;
  params.K = K;

  if nargin < 5
    params.varBoundBurnIn = 26;
    params.maxIteration     = 100;
  else
    params.varBoundBurnIn = burnin;
    params.maxIteration       = burnin * 2;
    if params.maxIteration < 100
      params.maxIteration  = 100;
    end
  end

  params.varBoundInterval = 5;
  params.varBoundConvergenceThresh  = 5e-4;

  params.lowerBound    = -inf;

  params.smoothKernelSize = 31;
  params.smoothKernelSigma = 11.0;
