function CBIG_KRDNN_KR_CV_example()

% CBIG_KRDNN_KR_CV_example()
% 
% This function is the example function of the kernel ridge regression
% algorithm with cross validation. 
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% add kernel regression code to path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath(genpath(fullfile(CBIG_CODE_DIR, '/utilities/matlab/predictive_models/KernelRidgeRegression/')));

% generate setup file
param = load(fullfile(CBIG_CODE_DIR, ...
	'stable_projects/preprocessing/Li2019_GSR/examples/output/KernelRidgeRegression/setup_file.mat'));
param.outdir = fullfile(pwd, 'output_kr_cv_example'); % Output folder of all output, can change to whatever you want
mkdir(param.outdir)
save(fullfile(param.outdir, 'setup_file.mat'), '-struct', 'param');

% run kernel regression
CBIG_KRR_workflow(fullfile(param.outdir, 'setup_file.mat'), 0);

% remove path
rmpath(genpath(fullfile(CBIG_CODE_DIR, '/utilities/matlab/predictive_models/KernelRidgeRegression/')));