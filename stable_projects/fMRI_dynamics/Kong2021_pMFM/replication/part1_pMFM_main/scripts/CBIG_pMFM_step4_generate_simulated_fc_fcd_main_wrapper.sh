#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step4_generate_simulated_fc_fcd.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../part1_pMFM_main/scripts
source activate pMFM
python CBIG_pMFM_step4_generate_simulated_fc_fcd.py

mv ../output/step4_MFM_simulated_data ../../replication/part1_pMFM_main/results/
