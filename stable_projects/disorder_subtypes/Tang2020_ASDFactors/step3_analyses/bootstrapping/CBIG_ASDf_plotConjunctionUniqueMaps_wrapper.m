function CBIG_ASDf_plotConjunctionUniqueMaps_wrapper(output_dir)
% CBIG_ASDf_plotConjunctionUniqueMaps_wrapper(output_dir)
% 
% Wrapper function to plot conjunction and unique maps for RSFC factors
%
% Input:
%     - output_dir:
%           Absolute path to directory where output results will be saved
%
% Example:
%       CBIG_ASDf_plotConjunctionUniqueMaps_wrapper('~/conjunction_uniq_maps')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
addpath(fullfile(CODE_DIR,'step3_analyses','bootstrapping'));

%% load pre-computed significant RSFC
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
INPUT_DIR = fullfile(UNIT_TEST_DIR,'results','bootstrapping');
load(fullfile(INPUT_DIR,'factor1_thresholded.mat'));
factor1 = corr_mat_masked;
load(fullfile(INPUT_DIR,'factor2_thresholded.mat'));
factor2 = corr_mat_masked;
load(fullfile(INPUT_DIR,'factor3_thresholded.mat'));
factor3 = corr_mat_masked;

%% binarize significant RSFC and sum across all factors
f1_bin = factor1 ~= 0;
f2_bin = factor2 ~= 0;
f3_bin = factor3 ~= 0;
counts = f1_bin + f2_bin + f3_bin;

%% unique map
uniq_map = counts;
uniq_map(uniq_map > 1) = 0;
uniq_map_f1 = factor1 .* uniq_map;
uniq_map_f2 = factor2 .* uniq_map;
uniq_map_f3 = factor3 .* uniq_map;
scalelim = [-1.6e-5, 1.6e-5];
CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(uniq_map_f1, scalelim, fullfile(output_dir,'uniq_map_F1'));
save(fullfile(output_dir,'uniq_map_F1.mat'), 'uniq_map_f1');
CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(uniq_map_f2, scalelim, fullfile(output_dir,'uniq_map_F2'));
save(fullfile(output_dir,'uniq_map_F2.mat'), 'uniq_map_f2');
CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(uniq_map_f3, scalelim, fullfile(output_dir,'uniq_map_F3'));
save(fullfile(output_dir,'uniq_map_F3.mat'), 'uniq_map_f3');

%% conjunction map, edges shared across 2 or 3 factors
conj_map = counts;
conj_map(conj_map < 2) = 0;
save(fullfile(output_dir,'conj_map.mat'), 'conj_map');

%% conjunction map, edges shared across all 3 factors
mask_conj = counts;
mask_conj(mask_conj < 3) = 0;
sum_abs = abs(factor1) + abs(factor2) + abs(factor3);
conj_map_all = mask_conj .* sum_abs;
scalelim = [-1e-4, 1e-4];
CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(conj_map_all, scalelim, ...
 fullfile(output_dir,'conj_map_allFactors'));
save(fullfile(output_dir,'conj_map_allFactors.mat'), 'conj_map_all');

%% Remove paths
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
rmpath(fullfile(CODE_DIR,'step3_analyses','bootstrapping'));
