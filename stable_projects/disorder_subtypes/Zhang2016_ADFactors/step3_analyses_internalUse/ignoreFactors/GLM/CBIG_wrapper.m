clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

%------- Which quantity? - Q
% 1 for MMSE, 2 for MEM, 3 for EF
Q = 3;

[betas, stats] = CBIG_fitGLM(SUBINFO_FILE, GMICV_FILE, Q);
baseline_acrossStages = CBIG_hypoTest_acrossStages(betas, stats);