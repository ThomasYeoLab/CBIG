#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step8_gene_expression_analysis_schaefer.m`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Schaefer100/scripts
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_pMFM_step8_gene_expression_analysis_schaefer;exit;"

mv ../output/CBIG_pMFM_step8_gene_expression_analysis \
../../../replication/part2_pMFM_control_analysis/Schaefer100_parcellation/results/
