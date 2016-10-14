function params = CBIG_olda_CreateEmptyParams(T, SEED, dict_size, corpus_size, sampling_functor, sampling_data, lambda_init)

% Create empty params struct of parameters for online LDA
% FORMAT params = CBIG_olda_CreateEmptyParams(T, SEED, dict_size, corpus_size, sampling_functor, sampling_data, lambda_init)
% 
% T                  = # of topics
% SEED               = seed used for Matlab's random number generator
% dict_size          = # of vocabulary words
% corpus_size        = # of documents
% sampling_functor   = sampling function
% sampling_data      = sampling data
% lambda_init        = initial value of lambda
%
% Written by B.T.Thomas Yeo and CBIG under MIT licence: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% random initialization
params.seed = SEED;
if(verLessThan('matlab', '7.6'))
    rand('twister', sum(10000*params.seed));
else
    RandStream.setDefaultStream(RandStream('mt19937ar','seed', 10000*params.seed));
end

params.T = T;
params.V = dict_size;
params.D = corpus_size;

%%%%%%%%%%%%%%%%%%%%%%%
% Outer EM convergence
%%%%%%%%%%%%%%%%%%%%%%%
params.max_iter             = 100;
params.pred_likelihood_iter = 200; % check every pred_likelihood_iter iterations


%%%%%%%%%%%%%%%%%%%%%%%%
% Inner EM convergence
%%%%%%%%%%%%%%%%%%%%%%%%
params.e_converge     = 1e-3;
params.max_e_iter     = 100;
params.m_converge     = 1e-5;
params.max_m_iter     = 1000;
params.estimate_hyper = 1; % 0 => fix alpha, eta; 1 =>  estimate


%%%%%%%%%%%%%%%%%%%%%%%
% Online parameters
%%%%%%%%%%%%%%%%%%%%%%%
params.batch_size = 100; % Note that hoffman says 500
params.tau = 1; % hoffman = 1
params.kappa = 0.9; % hoffman = 0.9; 0 means no online learning
params.samp_functor = sampling_functor;
params.samp_data = sampling_data;
if(params.kappa == 0) % turn off online learning
    params.batch_size = params.D;
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
params.init_alpha     = 50/T; % hyper on theta, same as ADU
params.init_eta       = 50/T; % hyper on beta, same as ADU
params.log_alpha      = log(params.init_alpha); 
params.log_eta        = log(params.init_eta); 

params.init_type = lambda_init;
if(strcmp(lambda_init, 'gamma'))
    params.lambda = gamrnd(100, 1/100, [params.T params.V]); % T x V
elseif(strcmp(lambda_init, 'rand'))
    params.lambda = rand([params.T params.V]); % T x V
elseif(strcmp(lambda_init, 'rand_docs'))
    doc_by_words = params.samp_functor(params, 'samp');
    params.lambda = doc_by_words(1:params.T, :) + rand([params.T params.V])/100;
else
    
end

% lambda ss
params.lambda_is_variational = 1;
if(params.lambda_is_variational)
    params.lambda_ss = bsxfun(@minus, digamma(params.lambda), digamma(sum(params.lambda, 2)));
else
    params.lambda = bsxfun(@times, params.lambda, 1./sum(params.lambda, 2)); % treat lambda as beta, normalize prob distribution
end



