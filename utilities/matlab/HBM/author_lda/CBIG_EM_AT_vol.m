function params = CBIG_EM_AT_vol(corpus, paradigm_by_exp, params, brain_mask)

% params = CBIG_EM_AT_vol(corpus, paradigm_by_exp, params, brain_mask)
%
% Runs EM algorithm to estimate parameters of the Author-Topic model
% Note that compared to CBIG_EM_AT.m, this function assumes the Author-Topic model is applied to activation
% foci in the brain volume
% FORMAT params   = CBIG_EM_AT_vol(corpus, paradigm_by_exp, params)
%
% corpus          = D x 1 cell where D is the number of documents
%                   corpus{d} is a Nd x V sparse matrix, where Nd is the number of activation foci (unique words)
%                   in the experiemnt (document), and V is the number of voxels. corpus{d}(n, j) = 1
% paradigm_by_exp = A x D logical sparse matrix, where A is the number of paradigms (authors),
%                   D is the number of experiments (documents).
%                   paradigm_by_exp(a, d) = 1 if paradigm (author) "a" is in experiment (document) d
% params          = struct of parameters used by the model
%                   params is created by CreateEmptyEM_ATParams
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

rand('twister', sum(10000*params.seed));

params = CBIG_EM_InitATParams_vol(corpus, paradigm_by_exp, params, brain_mask);
paradigm_by_exp = full(paradigm_by_exp);

params.log_likelihood = -inf;
tic;
for iter = 1:params.em_max_iter
    
    params.doc_log_likelihood = zeros(params.D, 1);
    params.new_theta = zeros([params.A params.T]);
    params.new_beta  = zeros([params.T params.V]);
    
    % e-step
    for d = 1:params.D
        if(size(corpus, 2) == 1)
            [params, params.doc_log_likelihood(d)] = CBIG_EM_doc_e_step(corpus{d}, paradigm_by_exp(:, d), params);
        elseif(size(corpus, 2) == 3)
            [params, params.doc_log_likelihood(d)] = CBIG_EM_doc_e_step_wc(corpus(d, :), paradigm_by_exp(:, d), params);
        else
            error('dimensions of corpus wrong');
        end
    end
    
    % m-step
    params.theta     = params.new_theta + params.virtual_theta;
    params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
    params.log_theta = log(params.theta);
    
    params.beta     = params.new_beta + params.virtual_beta;
    params.beta     = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
    params.log_beta = log(params.beta);

    if(max(params.theta(:)) >= params.max_theta)
	disp('Refactoring ...');
         
        paradigm_index = (max(params.theta, [], 2) >= params.max_theta);
        process_index  = (max(params.theta, [], 1) >= params.max_theta);

        params.theta(paradigm_index, :) = 1/params.T;
        params.theta(:, process_index)  = 1/params.T;         
    	params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
    	params.log_theta = log(params.theta);

        old_beta = params.beta;  
      	for t = 1:params.T
	    if(process_index(t))
	        params.beta(t, :) = sum(bsxfun(@times, old_beta(~process_index, :), rand([sum(~process_index) 1])), 1);	
	    end
        end
    	params.beta     = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
    	params.log_beta = log(params.beta);
        
        params.doc_log_likelihood(:) = -inf;
    end
    
    % compute convergence conditions
    new_log_likelihood = sum(params.doc_log_likelihood);
    if(isinf(params.log_likelihood))
        log_likelihood_improvement = inf;
    else
        log_likelihood_improvement = abs(new_log_likelihood - params.log_likelihood)/abs(params.log_likelihood);
    end
    params.log_likelihood = new_log_likelihood;
    disp(['Iter ' num2str(iter, '%03d') ': log likelihood = ' num2str(params.log_likelihood) ' (' num2str(log_likelihood_improvement) '), time elapsed = ' num2str(toc)]);
    tic;
    
    if(log_likelihood_improvement < params.em_convergence && iter >= params.em_min_iter && mod(iter, params.mod_smooth) ~= 1)
        break;
    end
    
    if(mod(iter, params.mod_smooth) == 0 && iter ~= params.em_max_iter)
        disp('Smoothing ...');
        
        % reshape
        process_maps = zeros([size(brain_mask.vol, 1)*size(brain_mask.vol, 2)*size(brain_mask.vol, 3) params.T]);
        for i = 1:params.T
            process_maps(brain_mask.vol(:) == 1, i) = params.beta(i, :);
        end
        process_maps = reshape(process_maps, [size(brain_mask.vol) params.T]);
        
        % smooth
        for i = 1:params.T
            for j = 1:params.init_smooth
                process_maps(:, :, :, i) = CBIG_Smooth3DVolumeWithMasks(squeeze(process_maps(:, :, :, i)), ...
                                                                   brain_mask.vol, 'SAME', 'box', 3);
            end
        end
        
        process_maps = reshape(process_maps, [size(process_maps, 1)*size(process_maps, 2)*size(process_maps, 3) size(process_maps, 4)]);
        process_maps = process_maps';
        process_maps = process_maps(:, brain_mask.vol(:) == 1);
        
        params.beta  = full(process_maps);
        params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
        params.log_beta = log(params.beta);
    end
end

% final e-step to compute final log likelihood
params.doc_log_likelihood = zeros(params.D, 1);
for d = 1:params.D
    if(size(corpus, 2) == 1)
        error('dimensions of corpus wrong');
    elseif(size(corpus, 2) == 3)
        q = CBIG_EM_doc_inference_wc(corpus(d, :), paradigm_by_exp(:, d), params);
        params.doc_log_likelihood(d) = CBIG_EM_doc_log_likelihood_wc(corpus(d, :), paradigm_by_exp(:, d), params, q);
    else
        error('dimensions of corpus wrong');
    end
end
params.log_likelihood = sum(params.doc_log_likelihood);
disp(['Final log likelihood: ' num2str(params.log_likelihood)]);

params = rmfield(params, 'log_theta');
params = rmfield(params, 'log_beta');
params = rmfield(params, 'new_theta');
params = rmfield(params, 'new_beta');


% rearrange beta into a volume again!
process_maps = zeros([size(brain_mask.vol, 1)*size(brain_mask.vol, 2)*size(brain_mask.vol, 3) params.T]);
for i = 1:params.T
    process_maps(brain_mask.vol(:) == 1, i) = params.beta(i, :);
end
process_maps = reshape(process_maps, [size(brain_mask.vol) params.T]);

beta = brain_mask;
beta.nframes = params.T;
beta.vol = process_maps;
params.beta = beta;

