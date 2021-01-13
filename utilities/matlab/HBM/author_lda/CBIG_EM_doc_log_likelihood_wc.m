function doc_log_likelihood = CBIG_EM_doc_log_likelihood_wc(w, paradigm, params, q) 

% doc_log_likelihood = CBIG_EM_doc_log_likelihood_wc(w, paradigm, params, q) 
%
% Compute log likelihood of a single experiment (document) given the updated auxillary parameter q.
% Note that compared to CBIG_EM_doc_log_likelihood_wc.m, this function takes in an argument w of a different structure
% FORMAT doc_log_likelihood = CBIG_EM_doc_log_likelihood_wc(w, paradigm, params)
%
% w           = 1 x 3 cell array
% w{1}        = Nd x V sparse matrix where Nd is the number of activation foci (unique
%               words) in the experiemnt (document), and V is the number of voxels
%               (vocabulary words). w{1} is the same as w{2} but has counts
% w{2}        = Nd x V sparse matrix where w{2}(n, v) = 1 if the n-th activation foci in i
%               the experiment is the v-th voxel
% w{3}        = 1 x Nd vector, where w{n} is the number of times the n-th unique activation
%               (word) in the experiment (document) appears
% Note that w{1} = bsxfun(@times, w{2}, w{3}');
% paradigm = A x 1 logical sparse matrix, where A is the number of paradigms (authors).
%            paradigm(a) = 1 if paradigm (author) "a" is in the current experiment (document)
% params   = struct of common parameters used by the model
% q        = Ad x T x Nd matrx, where Ad, and Nd are the number of paradigms (authors)
%            and unique activation foci (words) in the current experiment (document).
%            T is the number of components (topics)
%            q = p(author, topic, nth unique word | nth unique word, current estimate of beta and theta) 
%
% doc_log_likelihood = log likelihood of a single experiment (document)
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(w{3},1) ~= 1
    error('Input argument ''w{3}'' should be a row vector');
end
if size(paradigm,2) ~= 1
    error('Input argument ''paradigm'' should be a column vector');
end

log_likelihood_theta = sum(sum(sum(bsxfun(@times, bsxfun(@times, q, params.log_theta(paradigm, :)), ...
    reshape(w{3}, [1 1 length(w{3})])), 3), 2));

beta_update = squeeze(sum(q, 1)); % T x Nd*
if(length(w{3}) == 1)
    log_likelihood_beta  = sum(sum(bsxfun(@times, beta_update .* transpose(params.log_beta * w{2}'), ...
        reshape(w{3}, [1 length(w{3})])), 2));
else
    log_likelihood_beta  = sum(sum(bsxfun(@times, beta_update .* (params.log_beta * w{2}'), ...
        reshape(w{3}, [1 length(w{3})])), 2));
end

q_entropy = -sum(sum(sum(bsxfun(@times, q .* log(q), reshape(w{3}, [1 1 length(w{3})])), 3), 2));

doc_log_likelihood = log_likelihood_theta + log_likelihood_beta + q_entropy; 

