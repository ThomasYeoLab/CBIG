function params = InitAT_CVBParams(params, init_alpha, init_eta)
  % params = InitAT_CVBParams(params, init_alpha, init_eta)

  % params.Wb = cell array of size D x 1, where Wb{d} = Nd x V matrix
  % with Wb{d}(n, v) = 1 if the n-th unique word of document d is the
  % v-th vocabulary word
  %
  % params.Wc = cell array of size D x 1, where Wc{d} = Nd x 1 vector
  % with Wc{d}(n) = number of times the n-th unique word appears in the
  % d-th document
  %
  % doc_by_author = cell array of size D x 1, where doc_by_author{d} = Ad x A
  % matrix with doc_by_author{d}(l, a) = 1 if the l-th author in the d-th
  % document is author a

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  % load input data
  currDir                         = pwd;
  data                            = load(params.data_path);
  Wb                              = sparse(logical(data.act));   % D x V sparse logical matrix with Wb(d,v) = 1 if the v-th vocabulary  is in document d
  Wc                              = data.corpus; % Dx1 cell array with Wc{d}(n) = number of times the n-th unique word is in document d
  doc_by_author                   = sparse(logical(data.exp_by_paradigm));   % D x A sparse logical matrix with doc_by_author(d, a) = 1 if author a is in document d

  params.Wb                       = Wb;
  params.Wc                       = Wc;
  params.doc_by_author            = doc_by_author;
  params.doc_by_num_authors       = full(sum(doc_by_author,2));
  clear Wb Wc data;

  params.D                        = size(params.doc_by_author, 1);
  params.A                        = size(params.doc_by_author, 2);
  params.V                        = size(params.Wb, 2);

  % load brain mask
  brain_mask                      = MRIread(params.mask_path);
  if(sum(brain_mask.vol(:) == 1) ~= params.V)
     error('Brain mask is different size from dictionary'); 
  end
  clear brain_mask

  disp('Initializing parameters randomly');

  params.alpha                    = init_alpha;
  params.eta                      = init_eta;

  params.phi                      = cell(params.D, 1);

  params.total_words_count        = uint32(0);
  Nv                              = zeros([params.V 1], 'single'); % number of each vocabulary words in all documents
  params.Na                       = zeros([params.A 1], 'single'); % Na(a) = number of each vocabulary words in all documents that author a is inside
  
  params.author_indices           = false(params.D, params.A);
  params.word_indices             = cell(params.D, 1);
  params.word_counts              = zeros(params.D, 1, 'uint32');
  for d = 1:params.D
    Nd                          = size(params.Wc{d}, 1);
    Ad                          = params.doc_by_num_authors(d);
    params.phi{d}               = rand(Nd, params.K, Ad, 'single');
    params.phi{d}               = bsxfun(@times, params.phi{d}, 1./ sum(sum(params.phi{d},3), 2));

    params.total_words_count    = params.total_words_count + sum(params.Wc{d});
    logical_indices             = params.Wb(d,:);
    indices                     = find(logical_indices);
    Nv(indices)                 = Nv(indices) + params.Wc{d};
    for a = 1:params.A
      if doc_by_author(d, a) == 1
        params.Na(a)            = params.Na(a) + single(sum(params.Wc{d}));
      end;
    end;

    params.word_indices{d}      = logical_indices;
    params.word_counts(d)       = size(params.Wc{d}, 1);
    params.author_indices(d,:)  = params.doc_by_author(d,:);
  end;
  clear doc_by_author

  params.unique_Nv                = unique(Nv);
  params.word_indices_by_count    = cell(size(params.unique_Nv,1),1);
  for i = 1:size(params.unique_Nv, 1)
      params.word_indices_by_count{i} = find(Nv == params.unique_Nv(i));
  end;
  clear Nv;
end
