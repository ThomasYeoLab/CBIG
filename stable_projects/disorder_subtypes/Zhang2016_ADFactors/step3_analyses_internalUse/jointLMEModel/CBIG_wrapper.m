clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

K = 3;

%------- Which quantity? - Q
% 1 for MMSE, 2 for MEM, 3 for EF
Q = 3;

% Requirement for x-axis: auto, auto-center, or manual-center
xlimCmd = 'auto-center';

if K == 2
    stats = CBIG_fitLME_2t(SUBINFO_FILE, GMICV_FILE, Q);
    [~, decline] = hypoTest_2t(stats);
elseif K == 3
    stats = CBIG_fitLME_3t(SUBINFO_FILE, GMICV_FILE, Q);
    [~, decline] = hypoTest_3t(stats);
    decline_acrossStages = CBIG_hypoTest_3t_withinFactorAcrossStages(stats);
elseif K == 4
    stats = CBIG_fitLME_4t(SUBINFO_FILE, GMICV_FILE, Q);
    [~, decline] = hypoTest_4t(stats);
end

% Plot comparisons
system('rm *.eps');
% CBIG_plotCmp(baseline, xlimCmd);
CBIG_plotCmp(decline, xlimCmd);

decline_pvals = decline(:, 3);
save('decline_pvals.mat', 'decline_pvals');
