function [g, phi] = CBIG_olda_variational_inference(params, wc)

% Variational inference
% FORMAT [g, phi] = CBIG_olda_variational_inference(params, wc)
%
% wc     = 1 x V vector, where V is the vocabulary size, with wc(i) = # times word i appears
% g      = 1 x T vector gamma, where T is the number of topics.
% phi    = T x V
% lambda = T x V
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% compute document specific gamma and phi

if size(wc,1) ~= 1
    error('Input argument ''wc'' should be a row vector');
end

g = ones(1, params.T);
phi = zeros([params.T params.V]);
for i = 1:params.max_e_iter
    
    % update phi
    gamma_term       = digamma(g) - digamma(sum(g));
    if(params.lambda_is_variational)
        new_phi = exp(bsxfun(@plus, params.lambda_ss, gamma_term')); % treat lambda as variational parameter on beta
    else
        new_phi = bsxfun(@times, params.lambda, exp(gamma_term')); % treat lambda as beta
    end
    new_phi = bsxfun(@times, new_phi, 1./sum(new_phi, 1));
    
    % update gamma
    new_g = sum(bsxfun(@times, new_phi, wc), 2)' + exp(params.log_alpha);
    
    % change
    max_g_diff   = max(abs(g(:) - new_g(:)));
    max_phi_diff = max(abs(phi(:) - new_phi(:)));
    
    % update
    g = new_g;
    phi = new_phi;
    
    % break condition
    %disp(['Iter ' num2str(i, '%03d') ' : gamma change: ' num2str(max_g_diff) ', phi change: ' num2str(max_phi_diff)]);
    if(max_g_diff < params.e_converge && max_phi_diff < params.e_converge)
        break;
    end
end


if(isnan(sum(g)))
    error('gamma contains nan');
end


