function CBIG_ASDf_mean2mat(inputDir, outputDir, k, r, scalelim)
% CBIG_ASDf_mean2mat(inputDir, outputDir, k, r, scalelim)
%
% This function computes factor-specific RSFC patterns (i.e., E(RSFC patterns|Factor)), 
% plot it and save it as a 419x419 matrix.
%
% Input:
%     - inputDir:
%           Absolute path to the directory where the polarLDA estimation
%           results are saved
%     - outputDir:
%           Absolute path to the directory to save the 419x419 matrix
%           visualization plot
%     - k:
%           A string indicating the number of factors. e.g., k = '2'.
%     - r:
%           A string indicating the random initialization number. e.g., r = '30'.
%     - scalelim (optional):
%           Scale limit to plot the 419x419 matrix. If not specified,
%           scalelim will be set as [-maxAbs maxAbs] where maxAbs is the
%           max. absolute value in the matrix.
%
% Example:
%	CBIG_ASDf_mean2mat('~/example_output/estimate','~/example_output/visualizeFactors'
%   ,'3','20',[-1.6e-5 1.6e-5]) 
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%%Check input variables
if exist('scalelim','var')
    if size(scalelim, 1) > 1
        error('Input argument ''scalelim'' should be a row vector');
    end
end

inputDir = fullfile(inputDir, ['k' k], ['r' r]);
outputDir = fullfile(outputDir, ['k' k], ['r' r]);
num_ROIs = 419;

beta = load(fullfile(inputDir, 'final.beta'));
rho = load(fullfile(inputDir, 'final.rho'));

for factor_idx = 1:str2double(k)
    mean_corrmat = zeros(num_ROIs,num_ROIs);
    
    beta_row = exp(beta(factor_idx, :));
    rho_row = exp(rho(factor_idx, :));
    
    mean_row = beta_row.*(2*rho_row-1);
    
    % Get back to MxM matrix where M is the number of ROIs
    index = 0;
    for j = 1:(num_ROIs-1)
        for i = (j+1):num_ROIs
            index = index + 1;
            mean_corrmat(i,j) = mean_row(index);
        end
        mean_corrmat(i,i) = 0; % Diagonal set to 0
    end
    for i = 1:num_ROIs
        for j = (i+1):num_ROIs
            mean_corrmat(i,j) = mean_corrmat(j,i);
        end
    end
    
    save(fullfile(outputDir,['mean' num2str(factor_idx) '.mat']), 'mean_corrmat');
    
    % Plot the matrix
    if nargin <= 4 || isempty(scalelim)
        scalelim = [-max(abs(mean_corrmat(:))) max(abs(mean_corrmat(:)))];
    end
    CBIG_ASDf_Plot400Schaefer19Subcor17Networks_419by419Input(mean_corrmat, ...
        scalelim, fullfile(outputDir, ['mean' num2str(factor_idx)]));
    
end

