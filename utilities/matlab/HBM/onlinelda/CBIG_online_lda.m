function params = CBIG_online_lda(params)

% Implementation of online LDA
% FORMAT params = CBIG_online_lda(corpus, paradigm_by_exp, params)
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

e_beta = bsxfun(@times, params.lambda, 1./sum(params.lambda, 2));
last_iter = 0;
for iter = 1:params.max_iter
    
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute learning rate and sample data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.rho = (params.tau + iter - 1)^(-params.kappa);
    if(params.kappa == 0) % no online learning
        doc_by_words = params.samp_functor(params, 'full');
    else
        doc_by_words = params.samp_functor(params, 'samp');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E-step: update document specific gamma, phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.new_lambda = zeros([params.T params.V]);
    if(params.estimate_hyper)
        params.alpha_ss = 0;
    end
    
    for d = 1:size(doc_by_words, 1)
        params = CBIG_olda_e_step(params, doc_by_words(d, :));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E-step: update lambda (T x V)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    new_lambda = exp(params.log_eta) + (params.D/params.batch_size)*params.new_lambda;
    new_lambda = (1 - params.rho)*params.lambda + params.rho*new_lambda;
    if(params.lambda_is_variational)
        params.lambda_ss = bsxfun(@minus, digamma(new_lambda), digamma(sum(new_lambda, 2)));
        new_e_beta = bsxfun(@times, new_lambda, 1./sum(new_lambda, 2)); % compute expected beta
    else
        new_lambda = bsxfun(@times, new_lambda, 1./sum(new_lambda, 2)); % treat lambda as beta, so normalize as distribution
        new_e_beta = new_lambda;
    end
    lambda_change      = max(abs(new_lambda(:) - params.lambda(:)));
    lambda_frac_change = max(abs(new_lambda(:) - params.lambda(:))./params.lambda(:));
    beta_change        = max(abs(new_e_beta(:) - e_beta(:)));
    e_beta             = new_e_beta;
    params.lambda = new_lambda;
    disp_str = ['Iter ' num2str(iter, '%03d') ': lambda (' num2str(lambda_change) ', ' num2str(lambda_frac_change) ', ' num2str(beta_change) ')'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M-step for alpha and eta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(params.estimate_hyper)
        
        [params, log_alpha_change] = CBIG_olda_update_alpha(params);
        params.old_log_alpha(iter) = params.log_alpha;
        disp_str = [disp_str ', alpha = ' num2str(params.log_alpha) ' (' num2str(log_alpha_change) ')'];
        
        % only need to update eta if there is a prior on beta
        if(params.lambda_is_variational)
            [params, log_eta_change] = CBIG_olda_update_eta(params);
            params.old_log_eta(iter) = params.log_eta;
            disp_str = [disp_str ', eta = ' num2str(params.log_eta) ' (' num2str(log_eta_change) ')'];
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check log predictive prob
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(mod(iter, params.pred_likelihood_iter) == 0)
        
        log_predictive_likelihood = CBIG_olda_ComputeLogPredictiveLikelihood(params);
        
        if(last_iter == 0) % never been calculated before
            change_log_predictive_likelihood = inf;
        else
            change_log_predictive_likelihood = (log_predictive_likelihood - params.log_predictive_likelihood(end))/abs(params.log_predictive_likelihood(end));
        end
        last_iter = iter;
        params.log_predictive_likelihood(iter) = log_predictive_likelihood;
        disp_str = [disp_str ', likelihood = ' num2str(params.log_predictive_likelihood(iter)) ' (' num2str(change_log_predictive_likelihood) ')'];
        
        %save(['outputs/topic' num2str(params.T, '%03d') '/seed' num2str(params.seed, '%03d') '.mat'], 'params');
    end
    disp(disp_str);
    
    toc
end


if(iter ~= last_iter) % log_likelihood not computed in last iteration
    tic
    params.final_log_likelihood = CBIG_olda_ComputeLogPredictiveLikelihood(params);
    toc
else
    params.final_log_likelihood = params.log_predictive_likelihood(end);
end
disp(['Final Log Likelihood: ' num2str(params.final_log_likelihood)]);




if(0)
    tmp = bsxfun(@times, params.lambda, 1./sum(params.lambda, 2));
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex');
    cortex_label = [lh_avg_mesh.MARS_label == 2 rh_avg_mesh.MARS_label == 2];
    process = zeros(params.T, 20484);
    process(:, cortex_label) = tmp;
    for i = 1:params.T
        CBIG_DrawSurfaceMaps(process(i, 1:10242), process(i, 10243:end), 'fsaverage5', 'inflated', 0, max(process(i, :)));
        pause
        close all;
    end
end

