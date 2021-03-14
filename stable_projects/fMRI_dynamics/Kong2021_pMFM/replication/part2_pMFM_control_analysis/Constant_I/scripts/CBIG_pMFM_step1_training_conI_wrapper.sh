#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step1_training_conI.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Constant_I/scripts
source activate pMFM
python CBIG_pMFM_step1_training_conI.py

mv ../output/step1_training_results ../../../replication/part2_pMFM_control_analysis/Constant_I/results/
