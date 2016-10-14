function log_predictive_likelihood = CBIG_olda_ComputeLogPredictiveLikelihoodSingleDoc(params, obs_wc, unobs_wc)

% Compute log predictive likelihood of a single document
% FORMAT log_predictive_likelihood = CBIG_olda_ComputeLogPredictiveLikelihoodSingleDoc(params, obs_wc, unobs_wc)
% 
% obs_wc             = observed words
% unobs_wc           = unobserved words
%
% log_predictive_likelihood = log predictive likelihood of a single document
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


% expectation of theta conditioned on q(theta | gamma, obs_wc)
g = CBIG_olda_variational_inference(params, obs_wc);
theta = g/sum(g);

% expectation of beta conditioned on q(beta | lambda)
if(params.lambda_is_variational)
    beta = bsxfun(@times, params.lambda, 1./sum(params.lambda, 2)); 
else
    beta = params.lambda;
end

% Pr(w | w_obs, params)
pr_w = theta * beta;

% log predictive_likelhood
log_predictive_likelihood = sum(log(pr_w).*unobs_wc);  
