function [params, log_alpha_change] = CBIG_olda_update_alpha(params)

% Update alpha hyperparameter
% FORMAT [params, log_alpha_change] = CBIG_olda_update_eta(params)
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

old_log_alpha = params.log_alpha;
if(params.kappa == 0) % normal LDA
    
    init_alpha = params.init_alpha;
    for i = 1:params.max_m_iter
        
        % Compute gradient, hessian and update
        [log_alpha_gradient, alpha_gradient] = compute_log_alpha_gradient(params);
        log_alpha_hessian  = compute_log_alpha_hessian(params, alpha_gradient);
        log_alpha = params.log_alpha - log_alpha_gradient/log_alpha_hessian;
        
        % update alpha
        log_alpha_change = max(abs(log_alpha - params.log_alpha));
        params.log_alpha = log_alpha;
        
        % check exit condition and out of range
        if(isinf(exp(params.log_alpha)) || isnan(exp(params.log_alpha)))
            params.log_alpha = log(init_alpha);
            init_alpha = init_alpha*10;
            disp(['WARNING: (M-step opt) Alpha out of range! Reinitialize log(alpha) to ' num2str(params.log_alpha)]);
        else
            if(log_alpha_change < params.m_converge)
                break; 
            else
                disp(['Iter ' num2str(i, '%04d') ': Alpha change: ' num2str(log_alpha_change)]);
            end
        end
    end
    disp(['Final alpha: ' num2str(exp(params.log_alpha))]);
else
    % Compute gradient, hessian and update
    [log_alpha_gradient, alpha_gradient] = compute_log_alpha_gradient(params);
    log_alpha_hessian  = compute_log_alpha_hessian(params, alpha_gradient);
    log_alpha = params.log_alpha - params.rho*log_alpha_gradient/log_alpha_hessian;
    
    % update alpha
    params.log_alpha = log_alpha;
end
log_alpha_change = abs(params.log_alpha - old_log_alpha)/abs(old_log_alpha);

% check out of range condition
if(isinf(exp(params.log_alpha)) || isnan(exp(params.log_alpha)))
    disp(['WARNING: (M-step end) Alpha out of range! Reinitialize to ' num2str(params.init_alpha)]);
    params.log_alpha = log(params.init_alpha);
end





function [log_alpha_gradient, alpha_gradient] = compute_log_alpha_gradient(params)

% Note that there should be a multiplication D/batch size, but we ignore this term because it will cancel out with
% corresponding term in hessian
alpha_gradient     = params.alpha_ss + ...
    params.T*params.batch_size*(digamma(params.T*exp(params.log_alpha)) - digamma(exp(params.log_alpha)));
log_alpha_gradient = alpha_gradient*exp(params.log_alpha);


function log_alpha_hessian = compute_log_alpha_hessian(params, alpha_gradient)

% Note that there should be a multiplication D/batch size, but we ignore this term because it will cancel out with
% corresponding term in gradient
alpha_hessian      = params.T*params.batch_size*(params.T*psi(1, params.T*exp(params.log_alpha)) - psi(1, exp(params.log_alpha)));
log_alpha_hessian  = alpha_hessian*(exp(params.log_alpha)^2) + alpha_gradient*exp(params.log_alpha);
