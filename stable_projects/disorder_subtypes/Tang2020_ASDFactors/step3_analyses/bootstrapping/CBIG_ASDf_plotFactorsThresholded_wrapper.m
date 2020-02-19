function CBIG_ASDf_plotFactorsThresholded_wrapper(output_dir)
% CBIG_ASDf_plotFactorsThresholded_wrapper(output_dir)
%
% Wrapper function to plot bootstrap thresholded factors
%
% Input:
%     - output_dir:
%           Absolute path to directory where output results will be saved
%
% Example:
%       CBIG_ASDf_plotFactorsThresholded_wrapper('~/bootstrapping/visualization')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Load inputs
% Final threshold from FDR
final_thresh = 0.0186;

% Number of factors
k = 3;

% Scale limit for plots
scalelim = [-1.6e-5, 1.6e-5];

% p-values from estimated z-scores from bootstrapping
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile('stable_projects','disorder_subtypes','Tang2020_ASDFactors');
INPUT_DIR = fullfile(UNIT_TEST_DIR,'results','bootstrapping');
input_data = load(fullfile(INPUT_DIR,'bootstrappedByNetworks_pVals.mat'));
p_vals_blk = input_data.p_vals_blk;

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step2_polarLDA'));
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

%% Original factor
inputDir_best = fullfile(CODE_DIR,'data_release','files_polarLDA','model_K3');
beta_best = exp(load(fullfile(inputDir_best,'final.beta')));
rho_best = exp(load(fullfile(inputDir_best,'final.rho')));
Mean_best = beta_best.*(2*rho_best-1);
Mean_best = Mean_best([3 2 1],:); % re-order the original factors

%% Plot bootstrap thresholded factors 
for factor_idx = 1:k
    p = p_vals_blk(factor_idx,:);
    blk_indices = p <= final_thresh;
    blk_indices = blk_indices';
    
    filename_prefix = fullfile(output_dir, ['factor' num2str(factor_idx) '_thresholded']);
    corr_mat_masked = CBIG_ASDf_Plot400Schaefer19Subcor17Networks_thresholded(Mean_best(factor_idx,:), ...
scalelim, filename_prefix, blk_indices);
    close all;
    save([filename_prefix '.mat'], 'corr_mat_masked');
end

%% Remove paths
rmpath(fullfile(CODE_DIR,'step2_polarLDA'));
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
