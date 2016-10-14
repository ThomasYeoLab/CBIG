close all;
clear;
clc;

addpath(genpath('../../functions/'));

%% Read in Mem and EF

subInfo = csvread('../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv');
rid_dx = subInfo(:, [1 8]);

rid_viscode_examdate_mem_ef = CBIG_get_mem_ef(rid_dx(:, 1), 0); % all data of 810 subjects

%% Compute normalization parameters with baseline Mem and EF of CN subjects

ridList = cell2mat(rid_viscode_examdate_mem_ef(:, 1));
indBlCN = ismember(ridList, rid_dx(rid_dx(:, 2)==1, 1))... RID who is CN
    & strcmp(rid_viscode_examdate_mem_ef(:, 2), 'bl'); % baseline

tmp = rid_viscode_examdate_mem_ef(indBlCN, :); % variable for sanity check

mem_blCN = cell2mat(rid_viscode_examdate_mem_ef(indBlCN, 4));
memMean = mean(mem_blCN);
memStd = std(mem_blCN);

ef_blCN = cell2mat(rid_viscode_examdate_mem_ef(indBlCN, 5));
efMean = mean(ef_blCN);
efStd = std(ef_blCN);

%% Normalize all future data

rid_viscode_examdate_mem_ef_normalized = rid_viscode_examdate_mem_ef;

% Memory
mem_all = cell2mat(rid_viscode_examdate_mem_ef(:, 4));
rid_viscode_examdate_mem_ef_normalized(:, 4) = num2cell((mem_all-memMean)/memStd);

% EF
ef_all = cell2mat(rid_viscode_examdate_mem_ef(:, 5));
rid_viscode_examdate_mem_ef_normalized(:, 5) = num2cell((ef_all-efMean)/efStd);

save('normalizedMemEF.mat', 'rid_viscode_examdate_mem_ef_normalized');