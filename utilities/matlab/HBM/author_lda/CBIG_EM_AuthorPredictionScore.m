function [author_log_posterior_prob, author_classification_score] = CBIG_EM_AuthorPredictionScore(corpus, paradigm_by_exp, params)

% [author_log_posterior_prob, author_classification_score] = CBIG_EM_AuthorPredictionScore(corpus, paradigm_by_exp, params)
%
% Compute log posterior probability of the tasks (authors) assignment, and score of tasks (authors) classification given the estimated parameters.
% Note that as compared to CBIG_EM_AuthorPredictionScore_vol, this function takes in an argument corpus of a different structure
% FORMAT [author_log_posterior_prob, author_classification_score] = CBIG_EM_AuthorPredictionScore(corpus, paradigm_by_exp, params)
%
% corpus          = D x 1 cell where D is the number of documents
%                   corpus{d} is a Nd x V sparse matrix, where Nd is the number of activation foci (unique words)
%                   in the experiemnt (document), and V is the number of voxels. corpus{d}(n, j) = 1
% paradigm_by_exp = A x D logical sparse matrix, where A is the number of paradigms (authors),
%                   D is the number of experiments (documents).
%                   paradigm_by_exp(a, d) = 1 if paradigm (author) "a" is in experiment (document) d
% params          = struct of parameters used by the model
%
% author_log_posterior_prob = average log posterior probability of all experiments (documents)
% author_classification_score = average classification score of tasks (authors) across all documents
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

params.log_beta = log(params.beta);
params.log_theta = log(params.theta);
params.new_beta = zeros(size(params.beta));
params.new_theta = zeros(size(params.theta)); 

author_log_posterior_prob = zeros(length(corpus), 1);
author_classification_score = zeros(length(corpus), 1);
for d = 1:length(corpus)

    q = CBIG_EM_doc_inference(corpus{d}, logical(ones(size(paradigm_by_exp, 1), 1)), params);
    
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
