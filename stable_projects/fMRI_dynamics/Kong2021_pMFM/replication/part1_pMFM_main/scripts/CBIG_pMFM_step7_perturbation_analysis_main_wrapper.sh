#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step7_perturbation_analysis.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../part1_pMFM_main/scripts
source activate pMFM
python CBIG_pMFM_step7_perturbation_analysis.py

mv ../output/step7_perturbation_simulation ../../replication/part1_pMFM_main/results/
