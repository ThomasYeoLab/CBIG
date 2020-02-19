function [p, ADJ, NULL] = CBIG_ASDf_FCDiffAllSub_NBS(tThresh, id_asd_file, ...
id_con_file, sub_info_file, Nperm, output_name)
% [p, ADJ, NULL] = CBIG_ASDf_FCDiffAllSub_NBS(tThresh, id_asd_file, 
% id_con_file, sub_info_file, Nperm, output_name)
%
% This function compares FC difference between all ASD and all control participants. 
% Statistical significance is tested using network-based statistic (NBS; Zalesky et al., 2010). 
% The results will be saved in the directory specified by output_name.
% 
% Input:
%     - tThresh:
%           t-statistic threshold for NBS. E.g., 1.5 or 2
%     - id_asd_file:
%           Absolute path to the .mat file containing all ASD subjects'
%           IDs. Assuming the .mat file is named as "id_asd", and it is a
%           N1x1 cell array, where N1 is the number of ASD subjects, and the
%           i-th cell is the i-th ASD subject's ID (a string).
%     - id_con_file:
%           Absolute path to the .mat file containing all control subjects'
%           IDs. Assuming the .mat file is named as "id_con", and it is a
%           N2x1 cell array, where N2 is the number of control subjects, and
%           the i-th cell is the i-th control subject's ID (a string).
%     - sub_info_file:
%           Absolute path to the .csv file containing all subjects'
%           demographics information
%     - Nperm:
%           Number of permutations. E.g., 10000.
%     - output_name:
%           Output file name including the full path where the results of
%           NBS will be saved
% Output:
%     - p:
%           A vector of corrected p-values for each component identified in
%           NBS
%     - ADJ:
%           MxM matrix, where M is the number of edges in NBS input connectivity matrix
%           to NBS. Adjacency matrix identifying the edges comprising each
%           component. Edges corresponding to the first p-value stored in the 
%           vector p are assigned the value 1 in the adjacency matrix ADJ, 
%           edges corresponding to the second p-value are assigned the value 2, etc.
%     - NULL:
%           1xNperm matrix. A vector of Nperm samples from the null distribution of
%           maximal component size
%
% Example:
%     [p, ADJ, NULL] = CBIG_ASDf_FCDiffAllSub_NBS(2, id_asd, id_con,
%      '../examples/input/subInfo_inf.csv', 10000, 'nbs')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Read subjects' IDs
load(id_asd_file); % id_asd .mat file
load(id_con_file); % id_con .mat file

%% Generate permutation set using PALM package
permSet = CBIG_ASDf_genPermSetForNBS(id_asd, id_con, sub_info_file, Nperm);

fprintf('------Performing NBS for all ASD and control participants:\n');
fprintf('t-stats threshold:%f\n', tThresh);
fprintf('Number of permutations:%d\n', Nperm);

%% Load FC matrices
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
inputDir = fullfile(UNIT_TEST_DIR,'data');
load(fullfile(inputDir, 'corrMat_ASD_inf.mat')); % load ASD FC matrices
load(fullfile(inputDir, 'corrMat_Con_inf.mat')); % load control FC matrices

%% Run NBS
[p,ADJ,NULL] = CBIG_ASDf_NBS_permSetInput(corrMat_ASD, corrMat_Con, tThresh, permSet, 'both');

file_name = [output_name '.mat'];

save(file_name,'p','ADJ','NULL');

