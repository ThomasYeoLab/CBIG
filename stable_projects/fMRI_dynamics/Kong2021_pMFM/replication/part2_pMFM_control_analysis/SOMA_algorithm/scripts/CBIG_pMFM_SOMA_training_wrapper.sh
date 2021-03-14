#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_SOMA_training.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/SOMA_algorithm/scripts
source activate pMFM
python CBIG_pMFM_SOMA_training.py

mv ../output/low ../../../replication/part2_pMFM_control_analysis/SOMA_algorithm/results/
