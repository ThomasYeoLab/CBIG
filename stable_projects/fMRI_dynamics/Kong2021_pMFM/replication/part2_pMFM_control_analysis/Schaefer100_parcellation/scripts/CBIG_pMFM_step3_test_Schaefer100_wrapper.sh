#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step3_test_Schaefer100.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Schaefer100/scripts
source activate pMFM
python CBIG_pMFM_step3_test_Schaefer100.py

mv ../output/step3_test_results ../../../replication/part2_pMFM_control_analysis/Schaefer100_parcellation/results/
