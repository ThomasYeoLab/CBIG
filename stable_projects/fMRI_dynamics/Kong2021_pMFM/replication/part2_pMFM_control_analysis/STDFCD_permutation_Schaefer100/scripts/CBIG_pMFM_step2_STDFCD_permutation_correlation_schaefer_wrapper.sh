#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step2_STDFCD_permutation_correlation_schaefer.m`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/STDFCD_permutation_Schaefer100/scripts
matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all;CBIG_pMFM_step2_STDFCD_permutation_correlation_schaefer;exit;"

mv ../output/step2_compute_correlation \
../../../replication/part2_pMFM_control_analysis/STDFCD_permutation_Schaefer100/results/
