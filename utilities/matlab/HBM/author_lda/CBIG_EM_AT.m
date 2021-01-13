function params = CBIG_EM_AT(corpus, paradigm_by_exp, params)

% params = CBIG_EM_AT(corpus, paradigm_by_exp, params)
%
% Runs EM algorithm to estimate parameters of the Author-Topic model
% Note that compared to CBIG_EM_AT_vol.m, this function assumes the Author-Topic model is applied to activation
% foci on the brain surface
% FORMAT params   = CBIG_EM_AT(corpus, paradigm_by_exp, params)
%
% corpus          = D x 1 cell where D is the number of documents
%                   corpus{d} is a Nd x V sparse matrix, where Nd is the number of activation foci (unique words)
%                   in the experiemnt (document), and V is the number of voxels. corpus{d}(n, j) = 1
%                 or
% corpus          = 1 x 3 cell array
% corpus{1}       = Nd x V sparse matrix where Nd is the number of activation foci (unique
%                   words) in the experiemnt (document), and V is the number of voxels
%                  (vocabulary words). w{1} is the same as w{2} but has counts
% corpus{2}       = Nd x V sparse matrix where w{2}(n, v) = 1 if the n-th activation foci in i
%               the experiment is the v-th voxel
% corpus{3}       = 1 x Nd vector, where w{n} is the number of times the n-th unique activation
%               (word) in the experiment (document) appears
% Note that w{1} = bsxfun(@times, w{2}, w{3}');
%                   if nth word of document d is j-th vocabulary word.
% paradigm_by_exp = A x D logical sparse matrix, where A is the number of paradigms (authors),
%                   D is the number of experiments (documents).
%                   paradigm_by_exp(a, d) = 1 if paradigm (author) "a" is in experiment (document) d
% params          = struct of parameters used by the model
%                   params is created by CreateEmptyEM_ATParams
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


if size(corpus{3},1) ~= 1
    error('Input argument ''corpus{3}'' should be a row vector');
end


rand('twister', sum(10000*params.seed));
lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex');
cortex_label = [lh_avg_mesh.MARS_label == 2 rh_avg_mesh.MARS_label == 2];

params = CBIG_EM_InitATParams(corpus, paradigm_by_exp, params);
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
    
    % compute convergence conditions
    new_log_likelihood = sum(params.doc_log_likelihood);
    if(isinf(params.log_likelihood))
        log_likelihood_improvement = inf;
    else
        log_likelihood_improvement = abs(new_log_likelihood - params.log_likelihood)/abs(params.log_likelihood);
    end
    params.log_likelihood = new_log_likelihood;
    disp(['Iter ' num2str(iter, '%03d') ': log likelihood = ' num2str(params.log_likelihood) ...
        ' (' num2str(log_likelihood_improvement) '), time elapsed = ' num2str(toc)]);
    tic;
    
    if(log_likelihood_improvement < params.em_convergence && iter >= params.em_min_iter ...
            && mod(iter, params.mod_smooth) ~= 1)
        break;
    end
    
    if(mod(iter, params.mod_smooth) == 0 && iter ~= params.em_max_iter)
        disp('Smoothing ...');
        process_maps = zeros(params.T, length(lh_avg_mesh.MARS_label)+length(rh_avg_mesh.MARS_label));
        process_maps(:, cortex_label) = params.beta;
        process_maps(:, 1:10242) = MARS_AverageData(lh_avg_mesh, process_maps(:, 1:10242), 0, params.num_smooth);
        process_maps(:, 10243:end) = MARS_AverageData(lh_avg_mesh, process_maps(:, 10243:end), 0, params.num_smooth);
        
        process_maps = process_maps(:, cortex_label);
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
