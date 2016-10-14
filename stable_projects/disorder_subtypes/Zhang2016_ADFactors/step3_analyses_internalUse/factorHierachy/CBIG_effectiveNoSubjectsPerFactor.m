clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

%% Subjects

% All these info is baseline info
[RID_dx, ~, ~, ~, ~, ~] = get_data(SUBINFO_FILE, GMICV_FILE, 2);

% RIDs of each group
rid_AD = RID_dx(RID_dx(:, 2)==3, 1);
rid_ApMCI = CBIG_select_predementedGrp(RID_dx, 'a+_mci_147');
rid_ApCN = CBIG_select_predementedGrp(RID_dx, 'a+_hc_43');

%% K = 2

K = 2;

fprintf('\n\n=================== K = %d ===================\n\n', K);

% All these info is baseline info
[~, ~, rid_prob, ~, ~, ~] = get_data(SUBINFO_FILE, GMICV_FILE, K);

% What is the subtype order?
% Offset by 1, the RID column
TEM_COL = 1+2;
COR_COL = 1+1;
fprintf('\n\n\n');
fprintf('!!! Subtype order assumed to be: cortical, temporal+subcortical.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');

% AD
ind = ismember(rid_prob(:, 1), rid_AD);
fprintf('\n************** AD (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T+S:     %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in C:       %f\n', sum(rid_prob(ind, COR_COL)));

% A+ MCI
ind = ismember(rid_prob(:, 1), rid_ApMCI);
fprintf('\n************** A+ MCI (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T+S:     %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in C:       %f\n', sum(rid_prob(ind, COR_COL)));

% A+ CN
ind = ismember(rid_prob(:, 1), rid_ApCN);
fprintf('\n************** A+ CN (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T+S:     %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in C:       %f\n', sum(rid_prob(ind, COR_COL)));

%% K = 3

K = 3;

fprintf('\n\n=================== K = %d ===================\n\n', K);

% All these info is baseline info
[~, ~, rid_prob, ~, ~, ~] = get_data(SUBINFO_FILE, GMICV_FILE, K);

% What is the subtype order?
% Offset by 1, the RID column
TEM_COL = 1+2;
SUB_COL = 1+3;
COR_COL = 1+1;
fprintf('\n\n\n');
fprintf('!!! Subtype order assumed to be: cortical, temporal, subcortical.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');

% AD
ind = ismember(rid_prob(:, 1), rid_AD);
fprintf('\n************** AD (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T:       %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in S:       %f\n', sum(rid_prob(ind, SUB_COL)));
fprintf('Effective #subjects in C:       %f\n', sum(rid_prob(ind, COR_COL)));

% A+ MCI
ind = ismember(rid_prob(:, 1), rid_ApMCI);
fprintf('\n************** A+ MCI (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T:       %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in S:       %f\n', sum(rid_prob(ind, SUB_COL)));
fprintf('Effective #subjects in C:       %f\n', sum(rid_prob(ind, COR_COL)));

% A+ CN
ind = ismember(rid_prob(:, 1), rid_ApCN);
fprintf('\n************** A+ CN (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T:       %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in S:       %f\n', sum(rid_prob(ind, SUB_COL)));
fprintf('Effective #subjects in C:       %f\n', sum(rid_prob(ind, COR_COL)));

%% K = 4

K = 4;

fprintf('\n\n=================== K = %d ===================\n\n', K);

% All these info is baseline info
[~, ~, rid_prob, ~, ~, ~] = get_data(SUBINFO_FILE, GMICV_FILE, K);

% What is the subtype order?
% Offset by 1, the RID column
TEM_COL = 1+4;
SUB_COL = 1+2;
FRO_COL = 1+3;
PAR_COL = 1+1;
fprintf('\n\n\n');
fprintf('!!! Subtype order assumed to be: parietal, subcortical, frontal, temporal.\n');
fprintf('!!! Make sure this is true in the gamma file.\n');

% AD
ind = ismember(rid_prob(:, 1), rid_AD);
fprintf('\n************** AD (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T:       %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in S:       %f\n', sum(rid_prob(ind, SUB_COL)));
fprintf('Effective #subjects in F:       %f\n', sum(rid_prob(ind, FRO_COL)));
fprintf('Effective #subjects in P:       %f\n', sum(rid_prob(ind, PAR_COL)));

% A+ MCI
ind = ismember(rid_prob(:, 1), rid_ApMCI);
fprintf('\n************** A+ MCI (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T:       %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in S:       %f\n', sum(rid_prob(ind, SUB_COL)));
fprintf('Effective #subjects in F:       %f\n', sum(rid_prob(ind, FRO_COL)));
fprintf('Effective #subjects in P:       %f\n', sum(rid_prob(ind, PAR_COL)));

% A+ CN
ind = ismember(rid_prob(:, 1), rid_ApCN);
fprintf('\n************** A+ CN (N = %d) **************\n\n', sum(ind));
fprintf('Effective #subjects in T:       %f\n', sum(rid_prob(ind, TEM_COL)));
fprintf('Effective #subjects in S:       %f\n', sum(rid_prob(ind, SUB_COL)));
fprintf('Effective #subjects in F:       %f\n', sum(rid_prob(ind, FRO_COL)));
fprintf('Effective #subjects in P:       %f\n', sum(rid_prob(ind, PAR_COL)));