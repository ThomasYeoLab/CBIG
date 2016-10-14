clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

[rid_dx, ~, ~, ~, rid_sex, ~] = get_data(SUBINFO_FILE, GMICV_FILE, 3);

% Stats
maleAD = sum(rid_sex(rid_dx(:, 2)==3, 2)==1);
femaleAD = sum(rid_sex(rid_dx(:, 2)==3, 2)==2);
maleMCI = sum(rid_sex(rid_dx(:, 2)==2, 2)==1);
femaleMCI = sum(rid_sex(rid_dx(:, 2)==2, 2)==2);

% Table
%      Male    Female
% AD
% MCI
x = table([maleAD; maleMCI], [femaleAD; femaleMCI],...
    'VariableNames', {'Male', 'Female'}, 'RowNames', {'AD', 'MCI'});

% Fisher's exaact test
[~, p, stats] = fishertest(x); % introduced in R2014b
% p = 0.009