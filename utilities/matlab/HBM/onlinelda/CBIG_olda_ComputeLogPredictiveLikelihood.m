function log_predictive_likelihood = CBIG_olda_ComputeLogPredictiveLikelihood(params)

% Compute log predictive likelihood of all documents
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

log_predictive_likelihood = 0;
for i = 1:size(params.obs_wc, 1)
   log_predictive_likelihood = log_predictive_likelihood + ...
                               CBIG_olda_ComputeLogPredictiveLikelihoodSingleDoc(params, params.obs_wc(i, :), params.unobs_wc(i, :));
end

