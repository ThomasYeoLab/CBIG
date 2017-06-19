function [theta, beta] = EstimateATParams(params)
  % [theta, beta] = EstimateATParams(params)

  % Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

  theta   = ones(params.A, params.K, 'single') * params.alpha;
  beta    = ones(params.K, params.V, 'single') * params.eta;
  
  for d = 1:params.D
    word_indices    = params.word_indices{d};
    author_indices  = params.doc_by_author(d,:);
    Ad = params.doc_by_num_authors(d);
    
    tmp_1 = squeeze(sum(params.phi{d}, 1));
    if Ad > 1
        tmp_1 = tmp_1';   
    end;
    theta(author_indices, :) = theta(author_indices, :) + tmp_1;
    
    tmp_2 = squeeze(sum(params.phi{d}, 3));
    beta(:, word_indices) = beta(:, word_indices) + tmp_2';
  end;
  
  theta = bsxfun(@times, theta, 1./sum(theta, 2));
  beta = bsxfun(@times, beta, 1./sum(beta, 2));
