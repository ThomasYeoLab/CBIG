#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step6_SWSTD_state_main.m`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../part1_pMFM_main/scripts
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_pMFM_step6_SWSTD_state_main;exit;"

mv ../output/step6_SWSTD_state ../../replication/part1_pMFM_main/results/
