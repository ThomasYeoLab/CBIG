function CBIG_ASDf_visualizeFactors(inputDir, outputDir, k, scalelim)
% CBIG_ASDf_visualizeFactors(inputDir, outputDir, k, scalelim)
%
% Function to compute factor-specific RSFC patterns and visualize the
% factors for K = 2 and 3.
%
% Input:
%     - inputDir:
%           Absolute path to the directory where factor estimation results
%           are saved
%     - outputDir:
%           Absolute path to the directory where the factor-specific RSFC patterns 
%           (i.e., a 419x419 matrix for each factor) and the visualization 
%           plots will be saved.
%     - k:
%           String indicating the number of factors
%     - scalelim (optional):
%           Scale limit to plot the factor-specific RSFC patterns. If no
%           scalelim is entered, use the default [-1.6e-5 1.6e-5].
%
% Example:
%       CBIG_ASDf_visualizeFactors('~/Temporary/estimate',
%       '~/Temporary/visualizeFactors', '3', [-1.5e-5 1.5e-5]);
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%%Check input variables
if exist('scalelim','var')
    if size(scalelim, 1) > 1
        error('Input argument ''scalelim'' should be a row vector');
    end
end

%% Add path
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

%% Compute correlation between all the runs
% This function takes time if there are many runs
[corrRuns, ~, I] = CBIG_ASDf_computeCorrBetweenRuns(inputDir, k); 

%% Final estimate we obtained
r_final = num2str(I(1));
fprintf('Final estimate for K = %s: run %s.\n', k, r_final);
mkdir(fullfile(outputDir, ['k' k], ['r' r_final]));

% Plot correlation between final estimate and all other solutions
if size(corrRuns, 1) > 1
    CBIG_ASDf_plotCorrWithBest(inputDir, k, r_final, corrRuns, outputDir);
end

%% Plot the climbing of log-likelihood of the final estimate
% The log-likelihood should "go flat" if polarLDA has converged
CBIG_ASDf_plotLoglike(inputDir, k, r_final);


%% Compute E(RSFC patterns | Factor), project onto 419x419 matrix and plot it
if nargin <= 3 || isempty(scalelim)
    scalelim = [-1.6e-5 1.6e-5];
end
CBIG_ASDf_mean2mat(inputDir, outputDir, k, r_final, scalelim);

%% Write Pr(Factor | Participant), i.e., normalized \gamma in polarLDA, to .txt file, 
%% where each row is a participant's factor composition
gamma_file = fullfile(inputDir, ['k' k], ['r' r_final], 'final.gamma');
output_name = fullfile(outputDir, ['k' k], ['r' r_final], 'factorComp.txt');
CBIG_ASDf_gamma2table(gamma_file, output_name);

%% Remove path
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

end
