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
%  - params: updated struct with paramaters and auxillary variables.
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
  % Wb                              = sparse(logical(data.act));   % D x V sparse logical matrix with Wb(d,v) = 1 if the v-th vocabulary  is in document d
  % Wc                              = data.corpus; % Dx1 cell array with Wc{d}(n) = number of times the n-th unique word is in document d
  % doc_by_author                   = sparse(logical(data.exp_by_paradigm));   % D x A sparse logical matrix with doc_by_author(d, a) = 1 if author a is in document d
  Fb                              = sparse(logical(data.act));   % E x V sparse logical matrix with Fb(e,v) = 1 if voxel v is activated in experiment e
  Fc                              = data.corpus; % Ex1 cell array with Fc{e}(n) = number of times the n-th unique activation foci reported in experiment e

  % params.Wb                       = Wb;
  % params.Wc                       = Wc;
  % params.doc_by_author            = doc_by_author;
  % params.doc_by_num_authors       = full(sum(doc_by_author,2));
  % clear Wb Wc data;
  params.Fb                       = Fb;
  params.Fc                       = Fc;
  expByTask                     = sparse(logical(data.taskByExp'));
  params.expByTask              = expByTask;   % E x T sparse logical matrix with expByTask(e, t) = 1 if task t is utilized in experiment e
  params.expByNumTasks         = full(sum(params.expByTask, 2));
  clear Fb Fc data;

  % params.D                        = size(params.doc_by_author, 1);
  % params.A                        = size(params.doc_by_author, 2);
  % params.V                        = size(params.Wb, 2);
  params.E                        = size(params.expByTask, 1);
  params.T                        = size(params.expByTask, 2);
  params.V                        = size(params.Fb, 2);

  % load brain mask
  %   brain_mask                      = MRIread(params.mask_path);
  %   if(sum(brain_mask.vol(:) == 1) ~= params.V)
  %     error('Brain mask is different size from dictionary');
  %   end
  % clear brain_mask

  disp('Initializing parameters');

  params.alpha                    = alpha;
  params.eta                      = eta;

  params.phi                      = cell(params.E, 1);

  % params.total_words_count        = uint32(0);
  % Nv                              = zeros([params.V 1], 'single'); % number of each vocabulary words in all documents
  % params.Na                       = zeros([params.A 1], 'single'); % Na(a) = number of each vocabulary words in all documents that author a is inside
  params.totalFociCount         = uint32(0);
  Nv                              = zeros([params.V 1], 'single'); % number of times each brain voxel is activated across all experiments
  params.Nt                       = zeros([params.T 1], 'single'); % Nt(t) = number of activation times of all brain voxels reported in all experiments that utillize task t

  % params.author_indices           = false(params.D, params.A);
  % params.word_indices             = cell(params.D, 1);
  % params.word_counts              = zeros(params.D, 1, 'uint32');
  params.fociIndices             = cell(params.E, 1);
  params.fociCounts              = zeros(params.E, 1, 'uint32');
  for e = 1:params.E
    % Nd                          = size(params.Wc{d}, 1);
    % Ad                          = params.doc_by_num_authors(d);
    % params.phi{d}               = rand(Nd, params.K, Ad, 'single');
    % params.phi{d}               = bsxfun(@times, params.phi{d}, 1./ sum(sum(params.phi{d},3), 2));
    Ne                          = size(params.Fc{e}, 1);
    Te                          = params.expByNumTasks(e);
    params.phi{e}               = rand(Ne, params.K, Te, 'single');
    params.phi{e}               = bsxfun(@times, params.phi{e}, 1./ sum(sum(params.phi{e},3), 2));

    % params.total_words_count    = params.total_words_count + sum(params.Wc{d});
    % logicalIndices             = params.Wb(d,:);
    % indices                     = find(logicalIndices);
    % Nv(indices)                 = Nv(indices) + params.Wc{d};
    % for a = 1:params.A
    %   if doc_by_author(d, a) == 1
    %     params.Na(a)            = params.Na(a) + single(sum(params.Wc{d}));
    %   end;
    % end;
    params.totalFociCount     = params.totalFociCount + sum(params.Fc{e});
    logicalIndices             = params.Fb(e,:);
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
