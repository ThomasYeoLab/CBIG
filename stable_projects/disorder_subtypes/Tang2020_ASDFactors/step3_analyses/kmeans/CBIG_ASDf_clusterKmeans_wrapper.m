function CBIG_ASDf_clusterKmeans_wrapper(outputDir)
% CBIG_ASDf_clusterKmeans_wrapper(outputDir)
%
% Wrapper function to run k-means clustering on ASD participants' z-normalized RSFC.
%
% Input:
%     - outputDir:
%           Output directory to save the results.
%
% Example:
%     CBIG_ASDf_clusterKmeans('~/output/step1_FC2doc/zScores.mat', 3, 'correlation', 
%                               100, 1000, '~/output/kmeans')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Add paths
CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Tang2020_ASDFactors/step3_analyses/utilities']);
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Tang2020_ASDFactors/step2_polarLDA']);

%% Input data
UNIT_TEST_DIR = '/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/disorder_subtypes/Tang2020_ASDFactors';
INPUT_DIR = [UNIT_TEST_DIR '/results/results_long'];
zScores_dir = [INPUT_DIR '/FC2doc/step1_output_zScores.mat'];
k = 3;
distMet = 'correlation';
maxIter = 100;
replicates = 1000;

%% Set random seed
rng('default');

%% Run k-means
zScores = load(zScores_dir);
zScores = zScores.z;

[idx, C, sumD, D] = kmeans(zScores, k,'Distance', distMet, 'EmptyAction', 'singleton', ...
                            'MaxIter', maxIter, 'OnlinePhase', 'off', 'Replicates', replicates);
file_name = [outputDir '/kmeans_k=' num2str(k) '_distMet' distMet '_reps' num2str(replicates) ...
            '_maxIter' num2str(maxIter)];
save([file_name '.mat'], 'idx', 'C', 'sumD', 'D');

%% Plot cluster centroids
CBIG_ASDf_plotKmeans(C, file_name, [-1.8e-3 1.8e-3]);

%% Remove paths
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Tang2020_ASDFactors/step3_analyses/utilities']);
rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Tang2020_ASDFactors/step2_polarLDA']);
