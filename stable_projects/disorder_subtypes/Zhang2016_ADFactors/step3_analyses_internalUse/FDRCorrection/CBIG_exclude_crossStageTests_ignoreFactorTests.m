clear;
close all;
clc;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

% !!!!!!
% '_memRegEF' for regressing out each other; '' for not
tail = '';

%% MMSE baseline

cd('../jointGLM'); pwd

% K = 2
[betas, stats] = CBIG_fitGLM_2t(SUBINFO_FILE, GMICV_FILE, 1);
bl_mmse_k2 = CBIG_hypoTest_2t(betas, stats);
% K = 3
[betas, stats] = CBIG_fitGLM_3t(SUBINFO_FILE, GMICV_FILE, 1);
bl_mmse_k3 = CBIG_hypoTest_3t(betas, stats);
% K = 4
[betas, stats] = CBIG_fitGLM_4t(SUBINFO_FILE, GMICV_FILE, 1);
bl_mmse_k4 = CBIG_hypoTest_4t(betas, stats);

%% MMSE decline

cd('../jointLMEModel'); pwd

% K = 2
stats = CBIG_fitLME_2t(SUBINFO_FILE, GMICV_FILE, 1);
[~, decline_mmse_k2] = hypoTest_2t(stats);
% K = 3
stats = CBIG_fitLME_3t(SUBINFO_FILE, GMICV_FILE, 1);
[~, decline_mmse_k3] = hypoTest_3t(stats);
% K = 4
stats = CBIG_fitLME_4t(SUBINFO_FILE, GMICV_FILE, 1);
[~, decline_mmse_k4] = hypoTest_4t(stats);

%% MEM & EF baseline

cd(['../jointGLM' tail]); pwd

% MEM
% K = 2
[betas, stats] = CBIG_fitGLM_2t(SUBINFO_FILE, GMICV_FILE, 2);
bl_mem_k2 = CBIG_hypoTest_2t(betas, stats);
% K = 3
[betas, stats] = CBIG_fitGLM_3t(SUBINFO_FILE, GMICV_FILE, 2);
bl_mem_k3 = CBIG_hypoTest_3t(betas, stats);
% K = 4
[betas, stats] = CBIG_fitGLM_4t(SUBINFO_FILE, GMICV_FILE, 2);
bl_mem_k4 = CBIG_hypoTest_4t(betas, stats);

% EF
% K = 2
[betas, stats] = CBIG_fitGLM_2t(SUBINFO_FILE, GMICV_FILE, 3);
bl_ef_k2 = CBIG_hypoTest_2t(betas, stats);
% K = 3
[betas, stats] = CBIG_fitGLM_3t(SUBINFO_FILE, GMICV_FILE, 3);
bl_ef_k3 = CBIG_hypoTest_3t(betas, stats);
% K = 4
[betas, stats] = CBIG_fitGLM_4t(SUBINFO_FILE, GMICV_FILE, 3);
bl_ef_k4 = CBIG_hypoTest_4t(betas, stats);

%% MEM & EF decline

cd(['../jointLMEModel' tail]); pwd

% MEM
% K = 2
stats = CBIG_fitLME_2t(SUBINFO_FILE, GMICV_FILE, 2);
[~, decline_mem_k2] = hypoTest_2t(stats);
% K = 3
stats = CBIG_fitLME_3t(SUBINFO_FILE, GMICV_FILE, 2);
[~, decline_mem_k3] = hypoTest_3t(stats);
% K = 4;
stats = CBIG_fitLME_4t(SUBINFO_FILE, GMICV_FILE, 2);
[~, decline_mem_k4] = hypoTest_4t(stats);

% EF
% K = 2
stats = CBIG_fitLME_2t(SUBINFO_FILE, GMICV_FILE, 3);
[~, decline_ef_k2] = hypoTest_2t(stats);
% K = 3
stats = CBIG_fitLME_3t(SUBINFO_FILE, GMICV_FILE, 3);
[~, decline_ef_k3] = hypoTest_3t(stats);
% K = 4;
stats = CBIG_fitLME_4t(SUBINFO_FILE, GMICV_FILE, 3);
[~, decline_ef_k4] = hypoTest_4t(stats);

%% FDR over all tests

p = [bl_mmse_k2; bl_mmse_k3; bl_mmse_k4;...
    bl_mem_k2; bl_mem_k3; bl_mem_k4;...
    bl_ef_k2; bl_ef_k3; bl_ef_k4;...
    decline_mmse_k2; decline_mmse_k3; decline_mmse_k4;...
    decline_mem_k2; decline_mem_k3; decline_mem_k4;...
    decline_ef_k2; decline_ef_k3; decline_ef_k4];
p = p(:, 3);

alpha = 0.05
[~, finalThres] = FDR(p, alpha)

% '_memRegEF': 0.015555555555556