function CBIG_KRDNN_KR_TVT_example()

% CBIG_KRDNN_KR_TVT_example()
% 
% This function is the example function of the kernel ridge regression
% algorithm with normal training, validation and testing. 
% 
% Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

% add kernel regression code to path
addpath(genpath(fullfile(pwd, '/../KR_UKBB/')));

% run kernel regression
setup_file = '../examples/input/kr_tvt/mics_setup_example.txt';
lambda_set = [0.1 1 10];
CBIG_KRDNN_KRR_UKBB(setup_file, 4, lambda_set)

% remove path
rmpath(genpath(fullfile(pwd, '/../KR_UKBB/')));