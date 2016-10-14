clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

K = 3;

%------- Which quantity? - Q
% 1 for MMSE, 2 for MEM, 3 for EF
Q = 2;

% Requirement for x-axis: auto, auto-center, or manual-center
xlimCmd = 'auto-center';

if K == 2
    [betas, stats] = CBIG_fitGLM_2t(SUBINFO_FILE, GMICV_FILE, Q);
    baseline = CBIG_hypoTest_2t(betas, stats);
elseif K == 3
    [betas, stats] = CBIG_fitGLM_3t(SUBINFO_FILE, GMICV_FILE, Q);
    baseline = CBIG_hypoTest_3t(betas, stats);
    baseline_acrossStages = CBIG_hypoTest_3t_withinFactorAcrossStages(betas, stats);
elseif K == 4
    [betas, stats] = CBIG_fitGLM_4t(SUBINFO_FILE, GMICV_FILE, Q);
    baseline = CBIG_hypoTest_4t(betas, stats);
end

% Plot
system('rm *.eps');
CBIG_plotCmp(baseline, xlimCmd);

baseline_pvals = baseline(:, 3);
save('baseline_pvals.mat', 'baseline_pvals');

% Only for K = 3, generate an AD-only version of baseline for the paper
if K == 3
    baseline_AD = baseline(1:3, :);
    CBIG_plotCmp(baseline_AD, xlimCmd);
end