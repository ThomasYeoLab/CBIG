#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step6_SWSTD_state_Schaefer100.m`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Schaefer100/scripts
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_pMFM_step6_SWSTD_state_Schaefer100;exit;"

mv ../output/step6_SWSTD_state ../../../replication/part2_pMFM_control_analysis/Schaefer100_parcellation/results/
