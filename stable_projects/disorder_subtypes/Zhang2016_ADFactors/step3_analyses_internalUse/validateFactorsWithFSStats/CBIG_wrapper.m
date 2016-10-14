clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

% What is the subtype order?
TEM = 2;
SUB = 3;
COR = 1;
fprintf('\n\n\n');
fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');

% A structure is included if >= 50% its voxels have Pr(Voxel | Factor) > 0.8e-5
T = 0.5; 

% All these info is baseline info
[rid_dx, ~, rid_prob, ~, ~, ~] = get_data(SUBINFO_FILE, GMICV_FILE, 3);

% RID of 188 AD, 147 a+ MCI, 43 a+ AD
rid_188AD = rid_dx(rid_dx(:, 2)==3, 1);
rid_147MCI = CBIG_select_predementedGrp(rid_dx, 'a+_mci_147');
rid_43HC = CBIG_select_predementedGrp(rid_dx, 'a+_hc_43');
ridList = sort([rid_188AD; rid_147MCI; rid_43HC]);

%% Assign a factor to each structure

[label_structName_avgWinProb_tem, label_structName_avgWinProb_sub, label_structName_avgWinProb_cor] = ...
    CBIG_assignFactorsToStructures([TEM SUB COR]);

%% Sort included structures by average probabolity of the winning factor

structName_avgWinProb_tem = CBIG_sortByAvgProbWinningFactor(label_structName_avgWinProb_tem);
structName_avgWinProb_sub = CBIG_sortByAvgProbWinningFactor(label_structName_avgWinProb_sub);
structName_avgWinProb_cor = CBIG_sortByAvgProbWinningFactor(label_structName_avgWinProb_cor);
save('structName_avgWinProb.mat',...
    'structName_avgWinProb_tem', 'structName_avgWinProb_sub', 'structName_avgWinProb_cor');

%% Temporal

fprintf('\n************************ Temporal ************************\n\n');
rid_icv_vol = CBIG_getVol(ridList, label_structName_avgWinProb_tem(:, 2));
CBIG_fitGLM_hypoTest(rid_prob, rid_icv_vol, SUB+1, COR+1);

%% Subcortical

fprintf('\n************************ Subcortical ************************\n\n');
rid_icv_vol = CBIG_getVol(ridList, label_structName_avgWinProb_sub(:, 2));
CBIG_fitGLM_hypoTest(rid_prob, rid_icv_vol, SUB+1, COR+1);

%% Cortical

fprintf('\n************************ Cortical ************************\n\n');
rid_icv_vol = CBIG_getVol(ridList, label_structName_avgWinProb_cor(:, 2));
CBIG_fitGLM_hypoTest(rid_prob, rid_icv_vol, SUB+1, COR+1);