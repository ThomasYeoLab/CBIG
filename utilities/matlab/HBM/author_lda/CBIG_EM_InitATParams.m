function params = CBIG_EM_InitATParams(corpus, paradigm_by_exp, params)

% params = CBIG_EM_InitATParams(corpus, paradigm_by_exp, params)
%
% Initialize parameters for the EM algorithm of the Author-Topic model
% Note that compared to CBIG_EM_InitATParams_vol.m, this function allows the initialization of beta on the brain surface
% FORMAT params = CBIG_EM_InitATParams(corpus, paradigm_by_exp, params)
%
% corpus          = D x 1 cell where D is the number of documents
%                   corpus{d} is a Nd x V sparse matrix, where Nd is the number of activation foci (unique words)
%                   in the experiemnt (document), and V is the number of voxels. corpus{d}(n, j) = 1 
%                   if nth word of document d is j-th vocabulary word.
% paradigm_by_exp = A x D logical sparse matrix, where A is the number of paradigms (authors),
%                   D is the number of experiments (documents).
%                   paradigm_by_exp(a, d) = 1 if paradigm (author) "a" is in experiment (document) d
% params          = struct of parameters used by the model
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

params.A = size(paradigm_by_exp, 1);
params.V = size(corpus{1}, 2);
params.D = length(corpus);

if(strcmp(params.init, 'RAND'))
    disp('Initializing randomly');
    
    params.theta     = rand([params.A params.T]);
    params.theta     = params.theta;
    params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
    
    params.beta = rand([params.T params.V]);
    params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
elseif(strcmp(params.init, 'RAND_ACT'))
    
    params.theta     = rand([params.A params.T]);
    params.theta     = params.theta;
    params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
    
    act = zeros([params.A params.V]);
    for p = 1:params.A
        index = find(paradigm_by_exp(p, :) == 1 & sum(paradigm_by_exp, 1) == 1);
        
        for i = 1:length(index) % for each study in paradigm
            act(p, :) = act(p, :) + sum(corpus{index(i), 1}, 1);
        end
    end
    
    params.beta = zeros([params.T params.V]);
    for i = 1:params.T
        params.beta(i, :) = sum(bsxfun(@times, act, params.theta(:, i)), 1);
    end
    params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
else
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', 'fsaverage5', 'inflated', 'cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', 'fsaverage5', 'inflated', 'cortex');
    cortex_label = [lh_avg_mesh.MARS_label == 2 rh_avg_mesh.MARS_label == 2];

    if(strcmp(params.init, 'GIBBS'))
        
        disp('Initializing using Gibbs input');
        load(params.init_file);
        
        % initialize theta
        if(size(theta, 1) ~= params.A || size(theta, 2) ~= params.T)
            error('Initialization theta not the same size');
        end
        params.theta = full(theta);
        
        % initialize beta
        process_maps = process_maps';
        process_maps(:, 1:10242) = MARS_AverageData(lh_avg_mesh, process_maps(:, 1:10242), 0, params.init_smooth);
        process_maps(:, 10243:end) = MARS_AverageData(lh_avg_mesh, process_maps(:, 10243:end), 0, params.init_smooth);
        
        process_maps = process_maps(:, cortex_label);
        if(size(process_maps, 1) ~= params.T || size(process_maps, 2) ~= params.V)
            error('Initialization beta not the same size');
        end
        params.beta  = full(process_maps);
        params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
        
    elseif(strcmp(params.init, 'EM'))
        
        disp('Initializing using EM input');
        init_params = load(params.init_file, 'params');
        
        % initialize theta
        if(size(init_params.params.theta, 1) ~= params.A || size(init_params.params.theta, 2) ~= params.T)
            error('Initialization theta not the same size');
        end
        params.theta = init_params.params.theta;
        
        % initialize beta
        process_maps = zeros(params.T, length(lh_avg_mesh.MARS_label)+length(rh_avg_mesh.MARS_label));
        process_maps(:, cortex_label) = init_params.params.beta;
        process_maps(:, 1:10242) = MARS_AverageData(lh_avg_mesh, process_maps(:, 1:10242), 0, params.init_smooth);
        process_maps(:, 10243:end) = MARS_AverageData(lh_avg_mesh, process_maps(:, 10243:end), 0, params.init_smooth);
        
        process_maps = process_maps(:, cortex_label);
        if(size(process_maps, 1) ~= params.T || size(process_maps, 2) ~= params.V)
            error('Initialization beta not the same size');
        end
        params.beta  = full(process_maps);
        params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
        
    else
        error('Currently no other initialization');
    end
end

params.log_theta = log(params.theta); 
params.log_beta = log(params.beta);


