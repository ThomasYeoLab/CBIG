function [author_log_posterior_prob, author_classification_score] ...
    = CBIG_EM_AuthorPredictionScore_wc(corpus, paradigm_by_exp, params)

% [author_log_posterior_prob, author_classification_score] ...
% = CBIG_EM_AuthorPredictionScore_wc(corpus, paradigm_by_exp, params)
% Compute log posterior probability of the tasks (authors) assignment, and 
% score of tasks (authors) classification given the estimated parameters.
%
% Note that as compared to CBIG_EM_AuthorPredictionScore, this function 
% takes in an argument corpus of a different structure
% FORMAT [author_log_posterior_prob, author_classification_score] ...
% = CBIG_EM_AuthorPredictionScore_vol(corpus, paradigm_by_exp, params)
%
% corpus         = 1 x 3 cell array
% corpus{1}      = Nd x V sparse matrix where Nd is the number of activation foci (unique
%               words) in the experiemnt (document), and V is the number of voxels
%               (vocabulary words). w{1} is the same as w{2} but has counts
% corpus{2}      = Nd x V sparse matrix where w{2}(n, v) = 1 if the n-th activation foci in i
%               the experiment is the v-th voxel
% corpus{3}      = 1 x Nd vector, where w{n} is the number of times the n-th unique activation
%               (word) in the experiment (document) appears
% Note that corpus{1} = bsxfun(@times, corpus{2}, corpus{3}');
% paradigm_by_exp = A x D logical sparse matrix, where A is the number of paradigms (authors),
%                   D is the number of experiments (documents).
%                   paradigm_by_exp(a, d) = 1 if paradigm (author) "a" is in experiment (document) d
% params          = struct of parameters used by the model
%
% author_log_posterior_prob = average log posterior probability of all experiments (documents)
% author_classification_score = average classification score of tasks (authors) across all documents
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(corpus{3},1) ~= 1
    error('Input argument ''corpus{3}'' should be a row vector');
end

params.log_beta = log(params.beta);
params.log_theta = log(params.theta);

author_log_posterior_prob = zeros(length(corpus), 1);
author_classification_score = zeros(length(corpus), 1);
for d = 1:length(corpus)
    
    q = CBIG_EM_doc_inference_wc(corpus(d, :), logical(ones(size(paradigm_by_exp, 1), 1)), params);
    
    author_posterior_prob = mean(squeeze(sum(q, 2)), 2);
    
    author_log_posterior_prob(d) = mean(log(author_posterior_prob(paradigm_by_exp(:, d))));
    
    [Y, I] = sort(author_posterior_prob, 'descend');
    actual_paradigm = find(paradigm_by_exp(:, d));
    
    author_classification_score(d) = 0;
    for p = actual_paradigm'
        author_classification_score(d) = author_classification_score(d) + find(I == p);
    end
    author_classification(d) = author_classification_score(d)/length(actual_paradigm);
end

author_log_posterior_prob   = mean(author_log_posterior_prob);
author_classification_score = mean(author_classification_score);
