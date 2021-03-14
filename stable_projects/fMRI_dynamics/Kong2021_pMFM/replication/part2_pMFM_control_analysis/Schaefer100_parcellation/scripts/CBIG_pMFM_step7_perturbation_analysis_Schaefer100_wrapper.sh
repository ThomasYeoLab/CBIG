#!/bin/bash
# this function is the wrapper to run the `CBIG_pMFM_step7_perturbation_analysis.py`
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../../../../part2_pMFM_control_analysis/Schaefer100/scripts
source activate pMFM
python CBIG_pMFM_step7_perturbation_analysis.py

mv ../output/step7_perturbation_simulation \
../../../replication/part2_pMFM_control_analysis/Schaefer100_parcellation/results/
