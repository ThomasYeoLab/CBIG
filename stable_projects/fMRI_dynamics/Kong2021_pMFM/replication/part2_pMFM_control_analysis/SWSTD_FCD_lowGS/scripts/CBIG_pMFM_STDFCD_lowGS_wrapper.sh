#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_STDFCD_lowGS.m`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/SWSTD_FCD_lowGS/scripts
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_pMFM_STDFCD_lowGS();exit;"

mv ../output/step1_generate_order \
../../../replication/part2_pMFM_control_analysis/STDFCD_permutation_Schaefer100/results/
