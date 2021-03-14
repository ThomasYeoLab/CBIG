#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step1_training_consigma.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Constant_sigma/scripts
source activate pMFM
python CBIG_pMFM_step1_training_consigma.py

mv ../output/step1_training_results ../../../replication/part2_pMFM_control_analysis/Constant_sigma/results/
