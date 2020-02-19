function CBIG_ASDf_factorDistribution_wrapper(inputDir, outputDir)
% CBIG_ASDf_factorDistribution_wrapper(inputDir, outputDir)
%
% Wrapper function to visualize factor distribution for all ASD participants
% in ABIDE-II together, as well as in each site for both K = 2 and 3.
% 
% Input:
%     - inputDir:
%           Absolute path to the directory where the final estimate and factor visualization are located.
%           E.g., After you have run functions in step2_polarLDA and
%           obtained the latent factors, there will be a folder 'visualizeFactors' created. The full path to this
%           folder should be inputDir.
%     - outputDir:
%           Absolute path to the diectory where factor distribution plots will be saved
%
% Example:
%       CBIG_ASDf_factorDistribution_wrapper('~/example_output/visualizeFactors','~/example_output/analyses')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
CODE_DIR = fullfile(CBIG_CODE_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
addpath(fullfile(CODE_DIR,'step3_analyses','utilities'));

CBIG_REPDATA_DIR = getnev('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
DATA_DIR = fullfile(UNIT_TEST_DIR, 'data');
sub_info_file = fullfile(DATA_DIR, 'subInfo_654.csv');

%% plot factor distribution
k = '3';
r = '94';
factorOrder = [3 2 1];
factorComp = dlmread(fullfile(inputDir, ['k' k], ['r' r], 'factorComp.txt')); % K = 3 factor compositions
factorComp = factorComp(:,factorOrder);

figure;
CBIG_ASDf_visualizeFactorComp(sub_info_file, factorComp, str2double(k), fullfile(outputDir,'factorDistribution'));

%% Remove paths
rmpath(fullfile(CODE_DIR,'step3_analyses','utilities'));
