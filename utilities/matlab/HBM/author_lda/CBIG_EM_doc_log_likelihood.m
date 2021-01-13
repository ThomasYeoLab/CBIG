function doc_log_likelihood = CBIG_EM_doc_log_likelihood(w, paradigm, params, q) 

% doc_log_likelihood = CBIG_EM_doc_log_likelihood(w, paradigm, params, q) 
%
% Compute the log likelihood of a single experiment (document) given an updated
% auxillary parameter q.
% Note that compared to CBIG_EM_doc_log_likelihood_wc.m, this function takes in an argument w of a different structure
% FORMAT doc_log_likelihood = CBIG_EM_doc_log_likelihood(w, paradigm, params, q)
%
% w        = Nd x V sparse matrix, where Nd is the number of unique activation foci
%            (unique words) in the document and V is the number of voxels (vocabulary
%            words). w(n, v) = 1 if the n-th activatin foci (word) in the experiment
%            (document) is the v-th voxel (vocabulary word)
% paradigm = A x 1 logical sparse matrix, where A is the number of paradigms (authors).
%            paradigm(a) = 1 if paradigm (author) "a" is in the current experiment (document)
% params   = struct of common parameters used by the model
% q        = Ad x T x Nd matrx, where Ad, and Nd are the number of paradigms (authors)
%            and unique activation foci (words) in the current experiment (document).
%            T is the number of components (topics)
%            q = p(author, topic, nth word | nth word, current estimate of beta and theta) 
%
% doc_log_likelihood = log likelihood of a single experiment (document)
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(paradigm,2) ~= 1
    error('Input argument ''paradigm'' should be a column vector');
end

log_likelihood_theta = sum(sum(sum(bsxfun(@times, q, params.log_theta(paradigm, :)), 3), 2));

beta_update = squeeze(sum(q, 1)); % T x Nd
log_likelihood_beta  = sum(sum(beta_update .* (params.log_beta * w'), 2));

q_entropy = -sum(sum(sum(q .* log(q), 3), 2));

doc_log_likelihood = log_likelihood_theta + log_likelihood_beta + q_entropy; 
