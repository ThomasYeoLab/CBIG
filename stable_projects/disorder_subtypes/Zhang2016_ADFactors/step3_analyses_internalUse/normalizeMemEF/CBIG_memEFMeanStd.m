close all;
clear;
clc;

addpath(genpath('../../functions/'));

subInfo = csvread('../../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv');
rid_dx = subInfo(:, [1 8]);
ridList = rid_dx(:, 1);

rid_viscode_examdate_mem_ef = CBIG_get_mem_ef(ridList, 0);
indBl = strcmp(rid_viscode_examdate_mem_ef(:, 2), 'bl');
rid_mem_ef_bl = cell2mat(rid_viscode_examdate_mem_ef(indBl, [1 4 5]));

memBl = rid_mem_ef_bl(:, 2);
memBlMean = mean(memBl)
memBlStd = std(memBl)

efBl = rid_mem_ef_bl(:, 3);
efBlMean = mean(efBl)
efBlStd = std(efBl)