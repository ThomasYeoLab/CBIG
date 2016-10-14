clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

%!!!!!! Which score?
% 1 for MMSE, 2 for MEM, 3 for EF
Q = 2;

ps = CBIG_pairwiseCmpMeans(SUBINFO_FILE, GMICV_FILE, Q);