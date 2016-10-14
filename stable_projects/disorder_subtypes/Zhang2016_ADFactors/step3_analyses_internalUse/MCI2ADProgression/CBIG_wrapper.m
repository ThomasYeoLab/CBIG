close all;
clear;
clc;

addpath(genpath('../functions/'));

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

WINDOW_SIZE = 3; % years
K = 3;

system('rm *.eps');
CBIG_compare_progression(SUBINFO_FILE, GMICV_FILE, WINDOW_SIZE, K);
