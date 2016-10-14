clear;
clc;
close all;

addpath(genpath('../functions/'));

CBIG_plotSetup;

%% Assemble rid_prob and rid_ICVOverGM
r = dir('../../lda/postprocess/188blAD/k3/r*');

%----------------- bl
subInfo = csvread('../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv');
rid_gmVol_icv = csvread('../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv');
rid_ICVOverGM_810bl = [rid_gmVol_icv(:, 1), rid_gmVol_icv(:, 3)./rid_gmVol_icv(:, 2)];
% CN
rid_228CN = subInfo(subInfo(:, 8)==1, 1);
prob_228CN = CBIG_get_prob(['../../lda/outputs/k3' r.name '_inf228CN-gamma.dat']);
% MCI
rid_394MCI = subInfo(subInfo(:, 8)==2, 1);
prob_394MCI = CBIG_get_prob(['../../lda/outputs/k3' r.name '_inf394MCI-gamma.dat']);
% AD
rid_188AD = subInfo(subInfo(:, 8)==3, 1);
prob_188AD = CBIG_get_prob(['../../lda/postprocess/188blAD/k3/' r.name '/final.gamma']);
% Merge
rid_prob_810bl = [rid_228CN prob_228CN; rid_394MCI prob_394MCI; rid_188AD prob_188AD];

%----------------- m24
subInfo = csvread('../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_m24.csv');
rid_gmVol_icv = load('../../vbm_560m24/gmVolAndICV/rid_gmVol_icv_m24.csv');
rid_ICVOverGM_m24 = [rid_gmVol_icv(:, 1) rid_gmVol_icv(:, 3)./rid_gmVol_icv(:, 2)];
% Probs
rid_m24 = subInfo(:, 1);
prob_m24 = CBIG_get_prob(['../../lda/outputs/k3' r.name '_inf560m24-gamma.dat']);
assert(isequal(rid_m24,... % double check RID and probabilities are matched
    load('../../vbm_560m24/concat/gmToNonlinTemp_mod_4d_subjectOrder')));
rid_prob_m24 = [rid_m24 prob_m24];
% Ddivide by clinical group
rid_prob_m24_AD = rid_prob_m24(ismember(rid_prob_m24(:, 1), rid_188AD), :);
rid_prob_m24_MCI = rid_prob_m24(ismember(rid_prob_m24(:, 1), rid_394MCI), :);
rid_prob_m24_CN = rid_prob_m24(ismember(rid_prob_m24(:, 1), rid_228CN), :);

%% Test-retest

% Plot by cohort, color indicates amyloid status: red a+; green a-; blue unknown
% AD
rid_colors = CBIG_setColorByAmyloid(rid_prob_m24_AD(:, 1));
CBIG_testRetest_k3(rid_prob_810bl, rid_prob_m24_AD, rid_colors, 50, 'k');
% MCI
rid_colors = CBIG_setColorByAmyloid(rid_prob_m24_MCI(:, 1));
CBIG_testRetest_k3(rid_prob_810bl, rid_prob_m24_MCI, rid_colors, 50, 'k');
% CN
rid_colors = CBIG_setColorByAmyloid(rid_prob_m24_CN(:, 1));
CBIG_testRetest_k3(rid_prob_810bl, rid_prob_m24_CN, rid_colors, 50, 'k');

% Plot all together, color indicates clinical group: red AD; green CN; blue MCI
rid_prob_m24_all = [rid_prob_m24_AD; rid_prob_m24_MCI; rid_prob_m24_CN];
rid_colors = [num2cell(rid_prob_m24_all(:, 1)) repmat({'b'}, [size(rid_prob_m24_all, 1) 1])];
CBIG_testRetest_k3(rid_prob_810bl, rid_prob_m24_all, rid_colors, 25, 'r');