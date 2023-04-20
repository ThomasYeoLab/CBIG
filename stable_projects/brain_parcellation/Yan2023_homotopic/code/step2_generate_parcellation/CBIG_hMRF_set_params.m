function params = CBIG_hMRF_set_params(lhrh_avg_file, output_folder, varargin)
% params = CBIG_hMRF_set_params(lhrh_avg_file, output_folder, varargin);

% This function parses the input parameters into the output Matlab struct `params`.
% All input arguments are optional; if not specified, the argument would be configured according to their
% default values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-optional input parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - lhrh_avg_file: (character array)
%     Path to the premultiplied matrix. To generate the premultiplied matrices, refer to the functions under the 
%     folder ../step1_generate_fmri_input
%   - output_folder: (character array)
%     Path to the output folder.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configurable input parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - varargin: to pass in any parameter listed below with the format (paramName, paraVal). Refer to example usage
%     below. List of configurable parameters include:

%   Model Configuration:
%   - num_cluster: (integer) 
%     The total number of desirable number of clusters for the output parcellation. For instance, if the 
%     desired parcellation consists of 200 parcels on each side, num_cluster should be set to 400. We do not support
%     different numbers of parcels on the two hemispheres.
%   - mesh_type: (character array)
%     Currently supports 'fsaverage6'.
%   - grad_prior: (character array)
%     Currently only support 'gordon_fs6', which uses the gradient matrix derived in Generation and evaluation of 
%     a cortical area parcellation from resting-state correlations (Gordon et al., 2016)

%   Hyperparameters:
%   
%                                           *** Important note ***
%     Although default values are provided for the hyperparameters below, note that all parameters should also be tuned
%     based on the dimensionality of your input data (e.g., # of subjects, total # of time points.) The default values
%     below are good starting points if your input total time points is around 100k-500k (in our paper, the dimension of
%     the full GSP dataset is around 300k). In general, the smaller your input dimensionality is, the smaller your 
%     hyperparameters should be, and vice versa.
%
%   - initial_c: (double)
%     'initial_c' is the starting value of 'params.c'. 'params.c' will possibly change if `decrease_c` is turned on.
%     The hyperparameter c is for the gradient-weighted MRF term in the cost function. By default it is set to 80000.
%     The larger the value, the stronger the gradient-weighted MRF term would weigh in the cost function, thus
%     the parcellation tend to conform more to the gradient boundaries.
%   - initial_d: (double)
%     'initial_d' is the starting value of 'params.d'. 'params.d' will possibly change if `decrease_d` is turned on.
%     The hyperparameter that controls interhemispheric symmetry in the cost function. By default is set to 16000.
%     The stronger d is, the stronger the spatial symmetry constraint is (i.e., resultant parcellation more symmetric).
%   - k: (double)
%     The hyperparameter that controls the exponential decay in the gradient-weighted MRF term in the cost function.
%     By default you can set it to 15. The larger the value, the steeper the exponential decay. Refer to paper 
%     supplementary S2 for more mathematical details.
%   - w_xyz: (double)
%     The hyperparameter for the xyz terms that encourages parcel connectedness in the cost function. By default 
%     it is set to 100000. The stronger this term is, the rounder (i.e., stronger spatial connectedness) the parcels
%     would be (which may not be totally desirable). But if too weak, parcels tend to become disconnected.

%   Run Control:
%   - start_seed, num_rand_inits: (integer)
%     The program will generate parcellations based on random seeds [start_seed, start_seed + num_rand_inits).
%     With one core, one seed usually takes a few hours to run. It is recommended to set num_rand_inits to 1 
%     and run each seed in parallel.
%   - num_iterations: (integer)
%     The number of iterations to run the graphcut optimization algorithm with in 
%     CBIG_hMRF_update_labels_via_graphcut.
%   - premature_stopping: (logical)
%     The resultant parcellation from this program varies largely based on the random initialization. For certain 
%     random initializations, we could detect in earlier iterations that the resultant parcellation would not be
%     able to satisfy certain evaluation criteria. If `premature_stopping` is set to `true`, these "bad" random
%     initializations would be disgarded.
%     (fail fast principle).

%   Decrease-d algorithm:
%   - decrease_d: (logical)
%     If true, will run the decrease-d algorithm to relax the interhemispheric symmetry constraint.

%   Decrease-c algorithm:
%   - decrease_c: (logical)
%     If true, will run the decrease-c algorithm to relax the gradient-weighted MRF term.

%   Increase-tau algorithm setting:
%   - initial_tau: (double)
%     When set to zero, will initialize tau automatically via `CBIG_hMRF_initialize_lambda_in_vonmises_partition_func`. 
%     This is recommended since self-initialized value may cause numerical issues.
%   - increase_tau: (logical)
%     If true, will run increase-tau algorithm for parcels that are disconnected. When w_xyz is relatively large, 
%     this would not be an issue.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Non-configurable input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   - a: (integer)
%     The hyperparameter that controls interhemispheric symmetry in the d reduction step. Refer to paper 
%     supplementary for more details.
%   - nborhood_file: (character array)
%     Path to the matrix that defines the vertex-to-vertex neighborhood relationship between the two hemsipheres.
%     To replicate this neighborhood file, refer to ../utilities/CBIG_hMRF_generate_fs6_lhrh_nborhood.m
%   - lh_vert2rh_vert: (character array)
%     Path to the matrix that maps left hemisphere vertices to right hemisphere.
%   - rh_vert2lh_vert: (character array)
%     Path to the matrix that maps right hemisphere vertices to left hemisphere.
%   - CS: (character array)
%     Path to the file that stores the binary array indicating the location for central sulcus.
%   - V1: (character array)
%     Path to the file that stores the binary array indicating the location for V1.
%   - graphCutIterations: (integer)
%     The number of iterations within the graphcut optimization algorithm.
%   - increase_tau_iters: (integer)
%     Iterations by which we run the increase-tau algorithm. Alternatively, when all parcels become connected, we
%     would also terminate the increase-tau algorithm.
%   - tau_increase_speed: (integer)
%     The number of folds to increase tau by. Would recommend the default setting of 5.
%   - decrease_c_rate: (double)
%     The rate at which hyperparameter c is decreased within the decrease_c algorithm.
%   - min_verts_per_cluster: (integer)
%     Parcels containing equal or less than min_verts_per_cluster
%   - process_lost_parcels: (logical)
%     If true, when lost parcels are detected, they would be reassigned. For parcels that are lost on single
%     hemisphere, the singular parcel on the other hemi would be merged, before this parcel is reassigned on both 
%     hemispheres.
%   - flip_mismatched_parcels: (logical)
%     Mismatched parcels are those whose neighboring parcels are quite distinct on either hemispheres. If true, the
%     algorithm would try to flip pairs of mismatched parcels given that flipping would make both matched again.

% Output:
%   - params:
%     The output structure that specifies the input parameters.

% Example:
%   params = CBIG_hMRF_set_params(fullfile(inputdir, 'full_pmm.mat'),...
%     fullfile(output_dir, 'clustering'), 'start_seed', 1, 'num_rand_inits', 100,...
%    'num_cluster', 100, 'num_iterations', 100, 'initial_c', 1000, 'initial_d', 100, 'k', 15, 'w_xyz', 1000);
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

pnames = {'num_cluster' 'initial_c' 'initial_d' 'k' 'w_xyz' 'mesh_type' 'grad_prior' 'start_seed' 'num_rand_inits'...
 'num_iterations'  'premature_stopping' 'decrease_d'  'decrease_c'  'increase_tau'  'initial_tau'};
dflts = {400 80000 16000  15  100000  'fsaverage6'  'gordon_fs6'  1  1  100  true  true  true  false  0};

%% parse input arguments; if not given, use defaulVal.
[num_cluster, initial_c, initial_d, k, w_xyz, mesh_type, grad_prior, start_seed, num_rand_inits, num_iterations,...
 premature_stopping, decrease_d, decrease_c, increase_tau, initial_tau] = ...
 internal.stats.parseArgs(pnames, dflts, varargin{:});

%% configurable hyperparameters
% if passed in by bash, some input might be char instead of numerical. Need to make sure all input parameters are 
% correctly formatted.
getname = @(x) inputname(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-optional input parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(~exist('lhrh_avg_file', 'var'))
    error('Missing input parameter: lhrh_avg_file');
end
params.lhrh_avg_file = validationFcnChar(lhrh_avg_file, getname(lhrh_avg_file));
if(~exist('output_folder', 'var'))
    error('Missing input parameter: output_folder');
end
params.output_folder = validationFcnChar(output_folder, getname(output_folder));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configurable input parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.num_cluster = validationFcnInteger(num_cluster, getname(num_cluster));
params.initial_c = validationFcnNum(initial_c, getname(initial_c)); % save a copy of initial c
params.c = validationFcnNum(initial_c, getname(initial_c)); % this c would possibly be tuned
params.initial_d = validationFcnNum(initial_d, getname(initial_d)); % save a copy of initial d
params.d = validationFcnNum(initial_d, getname(initial_d)); % this d would possibly be tuned
params.k = validationFcnNum(k, getname(k));
params.w_xyz = validationFcnNum(w_xyz, getname(w_xyz));
params.mesh_type = validationFcnChar(mesh_type, getname(mesh_type));
params.grad_prior = validationFcnChar(grad_prior, getname(grad_prior));
params.start_seed = validationFcnInteger(start_seed, getname(start_seed));
params.num_rand_inits = validationFcnInteger(num_rand_inits, getname(num_rand_inits));
params.num_iterations = validationFcnInteger(num_iterations, getname(num_iterations));
params.premature_stopping = validationFcnLogical(premature_stopping, getname(premature_stopping));
params.decrease_d = validationFcnLogical(decrease_d, getname(decrease_d));
params.decrease_c = validationFcnLogical(decrease_c, getname(decrease_c));
params.increase_tau = validationFcnLogical(increase_tau, getname(increase_tau));
params.initial_tau = validationFcnNum(initial_tau, getname(initial_tau));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-configurable input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
input_data_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Yan2023_homotopic', 'code',...
 'step2_generate_parcellation', 'input');
params.a = 1;
params.nborhood_file = fullfile(input_data_dir, 'fsaverage6', 'derived_nborhood.mat');
params.lh_vert2rh_vert = fullfile(input_data_dir, 'fsaverage6', 'lh_to_rh.mat');
params.rh_vert2lh_vert = fullfile(input_data_dir, 'fsaverage6', 'rh_to_lh.mat');
params.CS = fullfile(input_data_dir, 'VGD11b_central_sulcus.mat');
params.V1 = fullfile(input_data_dir, 'V1_anchors.mat');
params.graphCutIterations = 1000;
params.increase_tau_iters = 1000;
params.tau_increase_speed = 5;
params.decrease_c_rate = 0.8;
params.min_verts_per_cluster = 3;
params.process_lost_parcels = true;
params.flip_mismatched_parcels = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(params.output_folder);

% configure convergence log file
params.convergence_log_folder = fullfile(params.output_folder, 'convergence_log');
mkdir(params.convergence_log_folder);

% Set gradient prior
if(strcmp(params.grad_prior, 'gordon_fs6'))
    params.lh_grad_file = fullfile(input_data_dir, 'fsaverage6', 'gordon_gradient_lh.mat');
    params.rh_grad_file = fullfile(input_data_dir, 'fsaverage6', 'gordon_gradient_rh.mat');
else
    error('Unknown gradient prior. Currently only support gordon gradient prior.\n');
end

% compute num of clusters PER HEMISPHERE
params.num_cluster_per_hemi = idivide(params.num_cluster, int32(2));

% format default output name, make output dir
params.output_name = [num2str(params.num_cluster),'parcels_C',...
    num2str(params.c, '%1.1e'), '_K', num2str(params.k), '_Wxyz',num2str(params.w_xyz,'%1.1e'), '_D',...
    num2str(params.d), '_A', num2str(params.a), '_iterations_', num2str(params.num_iterations)];

end

function output = validationFcnChar(input, varname)
% This function checks if input is char string
if(ischar(input))
    output = input;
else
    error('Input "%s" has to be a character string.', varname);
end
end

function output = validationFcnNum(input, varname)
% This function checks if input is numeric. If not, it will try to convert the input to number via str2num
% if str2num fails, an error would be thrown.
if(isnumeric(input))
    output = input;
else
    [out, ok] = str2num(input);
    if(ok)
        output = out;
    else
        error('Input "%s" has to be numeric or the input char string has to be convertible into a number.', varname);
    end
end
end

function output = validationFcnLogical(input, varname)
% This function checks if input is logical. If not, it will try to convert the input to logical via logical(str2num(x))
% if this fails, an error would be thrown.
if(islogical(input))
    output = input;
elseif(ischar(input))
    out = logical(str2num(input));
    if(isempty(out)) % if logical(x) gives empty, means conversion fails
        error('Input "%s" has to be logical, or the input char string should be convertible into a logical.', varname);
    else
        output = out;
    end
else
    error('Input "%s" has to be logical, or the input char string should be convertible into a logical.', varname);
end
end
    
function output = validationFcnInteger(input, varname)
% This function checks if input is integer, or a char string of an integer
input = validationFcnNum(input, varname); % first make sure input is a number / convertible to a number
if(rem(input, 1) ~= 0)
    error('Input "%s" has to be an integer', varname);
else
    output = input;
end
end
        