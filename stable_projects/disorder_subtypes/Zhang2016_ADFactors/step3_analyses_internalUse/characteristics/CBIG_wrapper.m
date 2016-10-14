clear;
clc;
close all;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

%!!!!!! 1 for a+ CN, 2 for a+ MCI, 3 for AD
dx = 1;

K = 3;

table = [];

% Baseline age
Q = 1;
table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];

if dx == 3
    % Onset age
    Q = 7;
    table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];
    
    % Year from onset to baseline
    Q = 6;
    table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];
end

% Edu
Q = 2;
table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];

% Gender
table = [table; CBIG_gender(SUBINFO_FILE, GMICV_FILE, K)];

% Amyloid
Q = 3;
table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];

% APOE-2
Q = 4;
table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];

% APOE-4
Q = 5;
table = [table; CBIG_blAge_edu_amyloid_apoe_disDura_onsetAge(SUBINFO_FILE, GMICV_FILE, K, Q, dx)];

save('table.mat', 'table');