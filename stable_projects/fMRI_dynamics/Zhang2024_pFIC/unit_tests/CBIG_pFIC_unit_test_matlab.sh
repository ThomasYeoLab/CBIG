#!/bin/bash
# This is a slightly modified version of the CBIG_pFIC_unit_test.sh to better suit the need of a MATLAB unit test
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

config=$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_config.ini
source activate pFIC

# training
python \
$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_training.py $config
python \
$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_transpose_matrix.py $config

conda deactivate

exit 0
