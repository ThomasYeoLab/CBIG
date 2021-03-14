#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step2_validation_T1T2.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/T1T2_only/scripts
source activate pMFM
python CBIG_pMFM_step2_validation_T1T2.py

mv ../output/step2_validation_results ../../../replication/part2_pMFM_control_analysis/T1T2_only/results/
