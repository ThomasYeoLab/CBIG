clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

K = 3;

% Requirement for x-axis: auto, auto-center, or manual-center
xlimCmd = 'auto-center';

if K == 2
    
elseif K == 3
    [betas, stats] = CBIG_fitGLM_3t(SUBINFO_FILE, GMICV_FILE);
    CBIG_hypoTest_3t(betas, stats);
elseif K == 4
    
end