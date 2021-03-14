#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_test_high_resolution.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/High_resolution/scripts
source activate pMFM
python CBIG_pMFM_test_high_resolution.py

mv ../output ../../../replication/part2_pMFM_control_analysis/High_resolution/results/
