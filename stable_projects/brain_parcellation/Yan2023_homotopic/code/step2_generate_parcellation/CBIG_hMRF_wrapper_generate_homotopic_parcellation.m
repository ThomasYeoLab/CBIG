function [results, params] = CBIG_hMRF_wrapper_generate_homotopic_parcellation(lhrh_avg_file, output_folder, varargin)
% [results, params] = CBIG_hMRF_wrapper_generate_homotopic_parcellation(lhrh_avg_file, output_folder, varargin)
%
% This is the wrapper function to set up and run the clustering algorithm based on the input arguments.

% For the notations below:
% N = no of vertices per hemisphere; 

% Input
%  - lhrh_avg_file: (character array)
%    Path to the 2Nx2N premultiplied matrix. To generate the premultiplied matrices, refer to the functions under the 
%    folder ./code/premultiplied_matrix
%  - output_folder: (character array)
%    Path to the output folder.

%  - varargin: (varargin)
%    Input format: ('variable_name1', variable_value1, 'variable_name2', variable_value2).
%    
%    You can input as many parameters as you want from the following list. If unspecified, default value would be used.
%    Each input variable must follow specific data type (e.g., Matlab logical, char array, double or integer).
%
%    For a starter, we recommend that the user should at least set the following parameters:
%    - num_cluster
%    - initial_c, initial_d, k, w_xyz (refer to our supplementary or replication wrappers for appropriate values)
%
%    You can also refer to CBIG_hMRF_set_params.m for more details on how the parameters are parsed.
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configurable input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%   - initial_c: (double)
%     'initial_c' is the starting value of 'params.c'. 'params.c' will possibly change if `decrease_c` is turned on.
%     The hyperparameter c is for the gradient-weighted MRF term in the cost function.
%     The larger the value, the stronger the gradient-weighted MRF term would weigh in the cost function, thus
%     the parcellation tend to conform more to the gradient boundaries.
%   - initial_d: (double)
%     'initial_d' is the starting value of 'params.d'. 'params.d' will possibly change if `decrease_d` is turned on.
%     The hyperparameter that controls interhemispheric symmetry in the cost function.
%     The stronger d is, the stronger the spatial symmetry constraint is (i.e., resultant parcellation more symmetric).
%   - k: (double)
%     The hyperparameter that controls the exponential decay in the gradient-weighted MRF term in the cost function.
%     The larger the value, the steeper the exponential decay. Refer to paper supplementary S2 for more details.
%   - w_xyz: (double)
%     The hyperparameter for the xyz terms that encourages parcel connectedness in the cost function.
%     The stronger this term is, the rounder (i.e., stronger spatial connectedness) the parcels
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
%      If true, will run the decrease-c algorithm to relax the gradient-weighted MRF term.

%   Increase-tau algorithm setting:
%   - initial_tau: (integer)
%     When set to zero, will use a customized function to initialize tau. This is recommended since self-initialized
%     value may cause numerical issues.
%   - increase_tau: (logical)
%     If true, will run increase-tau algorithm for parcels that are disconnected. When w_xyz is relatively large, 
%     this would not be an issue.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of configurable input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
%   - results: (struct)
%     A struct containing the final parcellation information, including most importantly:
%     + `lh_label`: the Nx1 full cortical label of the left hemisphere, with zeros indicating the medial wall vertices.
%     + `rh_label`: the Nx1 full cortical label of the right hemisphere, with zeros indicating the medial wall vertices.
%   - params: (struct)
%     A struct containing the parameters used in the parcellation generation algorithm.
%
% Example
%   - [results, params] = CBIG_hMRF_wrapper_generate_homotopic_parcellation(fullfile(inputdir, 'full_pmm.mat'),...
%     fullfile(output_dir, 'clustering'), 'start_seed', 1, 'num_rand_inits', 100,...
%    'num_cluster', 100, 'num_iterations', 100, 'initial_c', 1000, 'initial_d', 100, 'k', 15, 'w_xyz', 1000);
%
% Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation',...
 'Yan2023_homotopic', 'code', 'step2_generate_parcellation', 'lib'));

params = CBIG_hMRF_set_params(lhrh_avg_file, output_folder, varargin{:});
[params, results] = CBIG_hMRF_generate_parcellation_for_diff_rand_inits(params);

rmpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation',...
'Yan2023_homotopic', 'code', 'step2_generate_parcellation', 'lib'));
end
