function CBIG_ASDf_plotLoglike(inputDir, k, r)
% CBIG_ASDf_plotLoglike(inputDir, k, r)
%
% This function plots the "climbing" of the log-likelihood. If the model
% has converged, the plot should go flat.
% 
% Input:
%     - inputDir:
%           Absolute path to the input directory where the model estimation
%           results are saved
%     - k:
%           A string indicating the number of factors
%     - r:
%           A string indicating the run (random initialization) number
%
% Example:
%	CBIG_ASDf_plotLoglike('~/example_output/estimate','2','1')
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

figure;

log_likelihood = load(fullfile(inputDir, ['k' k], ['r' r], 'likelihood.dat'));

plot(log_likelihood(:, 1), 'o-');

xlabel('Iteration');
ylabel('Log-likelihood');

grid on;

