function params = CBIG_AuthorTopic_InitializeParams(params, alpha, eta)
% params = CBIG_AuthorTopic_InitializeParams(params, alpha, eta)
%
% Initialize the parameters of the Collapsed Variational Bayes (CVB)
% algorithm for the author-topic model
%
% Input:
%  - params     : struct containing parameters of the author-topic model
%                 and the CVB algorithm.
%  - alpha & eta: hyperparameters of the Dirichlet priors of the
%                 author-topic model
% Output:
%  - params: updated struct with parameters and auxiliary variables.
%
% Example:
%   params = CBIG_AuthorTopic_InitializeParams(params, 100, 0.01)
%   Initialize the parameters of the CVB algorithm with hyperparameters
%   alpha = 100 and eta = 0.01
%
% Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % load input data
  currDir                         = pwd;
  data                            = load(params.dataPath);
  Fb                              = sparse(logical(data.act));   % E x V sparse logical matrix with Fb(e,v) = 1 if voxel v is activated in experiment e
  Fc                              = data.corpus; % Ex1 cell array with Fc{e}(n) = number of times the n-th unique activation foci reported in experiment e

  params.Fb                       = Fb;
  params.Fc                       = Fc;
  expByTask                     = sparse(logical(data.taskByExp'));
  params.expByTask              = expByTask;   % E x T sparse logical matrix with expByTask(e, t) = 1 if task t is utilized in experiment e
  params.expByNumTasks         = full(sum(params.expByTask, 2));
  clear Fb Fc data;

  params.E                        = size(params.expByTask, 1);
  params.T                        = size(params.expByTask, 2);
  params.V                        = size(params.Fb, 2);


  disp('Initializing parameters');

  params.alpha                    = alpha;
  params.eta                      = eta;

  params.phi                      = cell(params.E, 1);

  params.totalFociCount         = uint32(0);
  Nv                              = zeros([params.V 1], 'single'); % number of times each brain voxel is activated across all experiments
  params.Nt                       = zeros([params.T 1], 'single'); % Nt(t) = number of activation times of all brain voxels reported in all experiments that utillize task t

  params.fociIndices             = cell(params.E, 1);
  params.fociCounts              = zeros(params.E, 1, 'uint32');
  for e = 1:params.E
    Ne                          = size(params.Fc{e}, 1);
    Te                          = params.expByNumTasks(e);
    params.phi{e}               = rand(Ne, params.K, Te, 'single');
    params.phi{e}               = bsxfun(@times, params.phi{e}, 1./ sum(sum(params.phi{e},3), 2));
    
    expFc                        = sum(params.Fc{e});
    assert(isfinite(expFc) && (floor(expFc) == expFc), ...
        ['Non-integer number of activate brain voxels  in exp ' num2str(e)])
    
    params.totalFociCount       = params.totalFociCount + expFc;
    logicalIndices              = params.Fb(e,:);
    indices                     = find(logicalIndices);
    Nv(indices)                 = Nv(indices) + params.Fc{e};
    for t = 1:params.T
      if expByTask(e, t) == 1
        params.Nt(t)            = params.Nt(t) + single(sum(params.Fc{e}));
      end;
    end;

    params.fociIndices{e}      = logicalIndices;
    params.fociCounts(e)       = size(params.Fc{e}, 1);
  end;
  clear expByTask

  params.uniqueNv                = unique(Nv);
  params.fociIndicesByCount    = cell(size(params.uniqueNv, 1), 1);
  for i = 1:size(params.uniqueNv, 1)
    params.fociIndicesByCount{i} = find(Nv == params.uniqueNv(i));
  end;
  clear Nv;
end
