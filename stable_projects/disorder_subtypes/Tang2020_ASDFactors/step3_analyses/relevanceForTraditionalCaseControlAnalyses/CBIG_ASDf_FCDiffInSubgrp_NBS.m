function [p, ADJ, NULL] = CBIG_ASDf_FCDiffInSubgrp_NBS(tThresh, id_asd_subgrp_file, ...
id_con_subgrp_file, id_asd_file, id_con_file, Nperm, output_name)
% [p, ADJ, NULL] = CBIG_ASDf_FCDiffInSubgrp_NBS(tThresh, id_asd_subgrp_file, 
% id_con_subgrp_file, id_asd_file, id_con_file, Nperm, output_name)
%
% This function compares FC difference between ASD and control participants
% in factor subgroups. Statistical significance is tested using network-based statistic 
% (NBS; Zalesky et al., 2010). 
% The results will be saved in the directory specified by output_name.
% 
% Input:
%     - tThresh:
%           t-statistic threshold for NBS. E.g., 1.5 or 2
%     - id_asd_subgrp_file:
%           Absolute path to the .mat file including all ASD subjects' IDs
%           in the subgroup. Assuming the .mat file is named as
%           "id_asd_subgrp", and it is a N_s1x1 cell array, where N_s1 is the
%           number of ASD subjects in the subgroup, and the i-th cell is
%           the i-th ASD subject's ID (a string)
%     - id_con_subgrp_file:
%           Absolute path to the .mat file including all control subjects' IDs
%           in the subgroup.. Assuming the .mat file is named as
%           "id_con_subgrp", and it is a N_s2x1 cell array, where N_s2 is the
%           number of control subjects in the subgroup, and the i-th cell is
%           the i-th control subject's ID (a string)
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
%     - Nperm:
%           Number of permutations. E.g., 10000.
%     - output_name:
%           Output file name including the full path where the results of
%           NBS will be saved
% Output:
%     - p:
%           1xC vector, where C is the number of identified components from 
%           NBS. Corrected p-values for each component from NBS
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
%     [p, ADJ, NULL] = CBIG_ASDf_FCDiffInSubgrp_NBS(1.5, id_asd_subgrp_file, 
%     id_con_subgrp_file, id_asd_file, id_con_file, 10000, 'nbs_subgroup')
%
% Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

%% Read subjects' IDs
load(id_asd_subgrp_file); % id_asd_subgrp .mat file
load(id_con_subgrp_file); % id_con_subgrp .mat file
load(id_asd_file); % id_asd .mat file
load(id_con_file); % id_con .mat file

%% Generate permutation set using PALM package
permSet = CBIG_ASDf_genPermSetForNBS(id_asd_subgrp, id_con_subgrp, Nperm);

fprintf('------Performing NBS in subgroups\n:');
fprintf('t-stats threshold:%f\n', tThresh);
fprintf('Number of permutations:%d\n', Nperm);

%% load FC matrices
CBIG_REPDATA_DIR = getenv('CBIG_REPDATA_DIR');
UNIT_TEST_DIR = fullfile(CBIG_REPDATA_DIR,'stable_projects','disorder_subtypes','Tang2020_ASDFactors');
inputDir = fullfile(inputDir, 'data');
load(fullfile(inputDir, 'corrMat_ASD_inf.mat')); % load ASD FC matrices
load(fullfile(inputDir, 'corrMat_Con_inf.mat')); % load control FC matrices

idx_asd_subgrp = ismember(id_asd, id_asd_subgrp);
idx_con_subgrp = ismember(id_con, id_con_subgrp);
corrMat_asd_subgrp = corrMat_ASD(:,:,idx_asd_subgrp);
corrMat_con_subgrp = corrMat_Con(:,:,idx_con_subgrp);

%% run NBS
[p,ADJ,NULL] = CBIG_ASDf_NBS_permSetInput(corrMat_asd_subgrp, corrMat_con_subgrp,tThresh, permSet, 'both');

file_name = [output_name '.mat'];

save(file_name,'p','ADJ','NULL');

