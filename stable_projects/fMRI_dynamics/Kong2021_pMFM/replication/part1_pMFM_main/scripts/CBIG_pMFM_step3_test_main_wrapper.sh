#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step3_test_main.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../part1_pMFM_main/scripts
source activate pMFM
python CBIG_pMFM_step3_test_main.py

mv ../output/step3_test_results ../../replication/part1_pMFM_main/results/
