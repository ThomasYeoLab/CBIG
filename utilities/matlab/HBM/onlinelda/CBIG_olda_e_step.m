function params = CBIG_olda_e_step(params, wc)

% E-step of the online LDA algorithm
% wc  = 1 x V vector, wc(i) = # times word i appear
% 
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if size(wc,1) ~= 1
    error('Input argument ''wc'' should be a row vector');
end

[g, phi] = CBIG_olda_variational_inference(params, wc);

% add document specific contribution to new lambda (lambda = T x V)
params.new_lambda = params.new_lambda + bsxfun(@times, phi, wc);

% add document specific contribution to alpha gradient
if(params.estimate_hyper)
    params.alpha_ss = params.alpha_ss + sum(digamma(g) - digamma(sum(g)));
end
