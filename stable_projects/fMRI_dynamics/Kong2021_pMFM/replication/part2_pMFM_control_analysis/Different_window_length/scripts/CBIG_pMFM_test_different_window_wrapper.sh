#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_test_different_window.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Different_window_length/scripts
source activate pMFM
python CBIG_pMFM_test_different_window.py

mv ../output/low ../../../replication/part2_pMFM_control_analysis/Different_window_length/results/
mv ../output/high ../../../replication/part2_pMFM_control_analysis/Different_window_length/results/
