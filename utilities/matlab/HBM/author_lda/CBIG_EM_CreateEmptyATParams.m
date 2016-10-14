function params = CBIG_EM_CreateEmptyATParams(SEED, T, init_type, init_file)

% params = CBIG_EM_CreateEmptyATParams(SEED, T, init_type, init_file)
%
% Generate a struct of parameters for the Author-Topic model
% FORMAT params = CBIG_EM_CreateEmptyATParams(SEED, T, init_type, init_file)
%
% SEED          = seed used by Matlab's random number generator
% T             = number of components (topics) to be estimated
% init_type     = a string used to determine how the auxillary paremeter
%                 of the model is generated. 'RAND' if the auxillary
%                 paramater is randomly generated, 'GIBBS' if the auxillary
%                 parameter is the output of Gibbs sampling
% init_file     = if init_type = 'GIBBS', the output of GIBBS sampling from
%                 init_file is used
%
% params        = struct of parameters used by the Author-Topic model
%
% Written by B.T.Thomas Yeo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

params.em_max_iter    = 1000;
params.em_min_iter    = 10;
params.em_convergence = 1e-5;
params.virtual_beta   = 0.01;
params.virtual_theta  = 50/T;
params.num_smooth     = 1;
params.mod_smooth     = 1000;
params.max_theta      = 1.1;

if(nargin < 3)
    params.init = 'RAND';
else
    params.init = init_type;
    params.init_file = init_file;
    params.init_smooth = 5;
end

params.seed = SEED;
params.T = T;
