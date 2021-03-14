#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step1_training_gradient.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Gradient_only/scripts
source activate pMFM
python CBIG_pMFM_step1_training_gradient.py

mv ../output/step1_training_results ../../../replication/part2_pMFM_control_analysis/Gradient_only/results/
