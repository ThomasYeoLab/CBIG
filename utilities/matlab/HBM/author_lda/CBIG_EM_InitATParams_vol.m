function params = CBIG_EM_InitATParams_vol(corpus, paradigm_by_exp, params, brain_mask)

% params = CBIG_EM_InitATParams_vol(corpus, paradigm_by_exp, params, brain_mask)
%
% Initialize parameters for the EM algorithm of the Author-Topic model
% Note that compared to CBIG_EM_InitATParams.m, this function initializes
% beta inside a brain volume
% FORMAT params = CBIG_EM_InitATParams_vol(corpus, paradigm_by_exp, params)
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

if(sum(brain_mask.vol(:) == 1) ~= params.V)
   error('Brain mask is different size from dictionary'); 
end

if(strcmp(params.init, 'RAND'))
    disp('Initializing randomly');
    
    params.theta     = rand([params.A params.T]);
    params.theta     = params.theta;
    params.theta     = bsxfun(@times, params.theta, 1./sum(params.theta, 2));
    
    params.beta = rand([params.T params.V]);
    params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
elseif(strcmp(params.init, 'RAND_ACT'))
    
    
else
    

    if(strcmp(params.init, 'GIBBS'))
        
        disp('Initializing using Gibbs input');
        load(params.init_file);
        
        % initialize theta
        if(size(theta, 1) ~= params.A || size(theta, 2) ~= params.T)
            error('Initialization theta not the same size');
        end
        params.theta = full(theta);
        
        % initialize beta
        for i = 1:params.T
            for j = 1:params.init_smooth
                process_maps.vol(:, :, :, i) = CBIG_Smooth3DVolumeWithMasks(squeeze(process_maps.vol(:, :, :, i)), ...
                                                                   brain_mask.vol, 'SAME', 'box', 3);
            end
        end
        
        process_maps = reshape(process_maps.vol, [size(process_maps.vol, 1)*size(process_maps.vol, 2)*size(process_maps.vol, 3) size(process_maps.vol, 4)]);
        process_maps = process_maps';
        process_maps = process_maps(:, brain_mask.vol(:) == 1);
        
        if(size(process_maps, 1) ~= params.T || size(process_maps, 2) ~= params.V)
            error('Initialization beta not the same size');
        end
        params.beta  = full(process_maps);
        params.beta = bsxfun(@times, params.beta, 1./sum(params.beta, 2));
        
    elseif(strcmp(params.init, 'EM'))
 
        
    else
        error('Currently no other initialization');
    end
end

params.log_theta = log(params.theta); 
params.log_beta = log(params.beta);


