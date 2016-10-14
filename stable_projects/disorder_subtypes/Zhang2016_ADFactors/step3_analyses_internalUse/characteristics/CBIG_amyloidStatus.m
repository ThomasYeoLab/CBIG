clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

[rid_dx, ~, ~, ~, ~, ~] = get_data(SUBINFO_FILE, GMICV_FILE, 3);

rid_AD = rid_dx(rid_dx(:, 2)==3, 1);
rid_MCI = rid_dx(rid_dx(:, 2)==2, 1);
rid_CN = rid_dx(rid_dx(:, 2)==1, 1);

%% AD
rid_dummy_dummy_amyloid = CBIG_get_amyloid(rid_AD);
rid_amyloid = cell2mat(rid_dummy_dummy_amyloid(:, [1 4]));

fprintf('\nAD total: %d\n', numel(rid_AD));
fprintf('AD having amyloid: %d\n', size(rid_amyloid, 1));
fprintf('AD amyloid+: %d\n', numel(rid_amyloid(rid_amyloid(:, 2)<192, 1)));

%% MCI
rid_dummy_dummy_amyloid = CBIG_get_amyloid(rid_MCI);
rid_amyloid = cell2mat(rid_dummy_dummy_amyloid(:, [1 4]));

fprintf('\nMCI total: %d\n', numel(rid_MCI));
fprintf('MCI having amyloid: %d\n', size(rid_amyloid, 1));
fprintf('MCI amyloid+: %d\n', numel(rid_amyloid(rid_amyloid(:, 2)<192, 1)));

%% AD
rid_dummy_dummy_amyloid = CBIG_get_amyloid(rid_CN);
rid_amyloid = cell2mat(rid_dummy_dummy_amyloid(:, [1 4]));

fprintf('\nCN total: %d\n', numel(rid_CN));
fprintf('CN having amyloid: %d\n', size(rid_amyloid, 1));
fprintf('CN amyloid+: %d\n', numel(rid_amyloid(rid_amyloid(:, 2)<192, 1)));