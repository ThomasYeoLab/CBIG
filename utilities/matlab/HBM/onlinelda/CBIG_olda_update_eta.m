function [params, log_eta_change] = CBIG_olda_update_eta(params)

% Update eta hyperparameter
% FORMAT [params, log_eta_change] = CBIG_olda_update_eta(params)
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

old_log_eta = params.log_eta;
if(params.kappa == 0) % normal LDA
    
    init_eta = params.init_eta;
    for i = 1:params.max_m_iter
        
        % Compute gradient, hessian and update
        [log_eta_gradient, eta_gradient] = compute_log_eta_gradient(params);
        log_eta_hessian  = compute_log_eta_hessian(params, eta_gradient);
        log_eta = params.log_eta - log_eta_gradient/log_eta_hessian;
        
        % update eta
        log_eta_change = max(abs(log_eta - params.log_eta));
        params.log_eta = log_eta;
        
        % check exit condition and out of range
        if(isinf(exp(params.log_eta)) || isnan(exp(params.log_eta)))
            params.log_eta = log(init_eta);
            init_eta = init_eta*10;
            disp(['WARNING: (M-step opt) Eta out of range! Reinitialize log(eta) to ' num2str(params.log_eta)]);
        else
            if(log_eta_change < params.m_converge)
                break;
            else
                disp(['Iter ' num2str(i, '%04d') ': Eta change: ' num2str(log_eta_change)]);
            end
        end
    end
else
    % Compute gradient, hessian and update
    [log_eta_gradient, eta_gradient] = compute_log_eta_gradient(params);
    log_eta_hessian = compute_log_eta_hessian(params, eta_gradient);
    log_eta = params.log_eta - params.rho*log_eta_gradient/log_eta_hessian;
    
    % update eta
    params.log_eta = log_eta;
end
log_eta_change = abs(params.log_eta - old_log_eta)/abs(old_log_eta);


% check out of range condition
if(isinf(exp(params.log_eta)) || isnan(exp(params.log_eta)))
    disp(['WARNING: (M-step end) Eta out of range! Reinitialize to ' num2str(params.init_eta)]);
    params.log_eta = log(params.init_eta);
end





function [log_eta_gradient, eta_gradient] = compute_log_eta_gradient(params)

eta_gradient     = sum(sum(bsxfun(@minus, digamma(params.lambda), digamma(sum(params.lambda, 2))))) + ...
    params.T*params.V*(digamma(params.V*exp(params.log_eta)) - digamma(exp(params.log_eta)));
log_eta_gradient = eta_gradient*exp(params.log_eta);


function log_eta_hessian = compute_log_eta_hessian(params, eta_gradient)

eta_hessian     = params.T*params.V*(params.V*psi(1, params.V*exp(params.log_eta)) - psi(1, exp(params.log_eta)));
log_eta_hessian = eta_hessian*(exp(params.log_eta)^2) + eta_gradient*exp(params.log_eta);
