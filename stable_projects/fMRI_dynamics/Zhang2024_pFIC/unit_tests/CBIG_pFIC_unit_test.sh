#!/bin/bash
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

MATLAB=`which matlab`

config=$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_config.ini
source activate pFIC

# training
python \
$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_training.py $config
python \
$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_transpose_matrix.py $config

conda deactivate

# check results against reference output
util_path=$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/examples
user_output_path=$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/examples/output
reference_output_path=$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/examples/reference_output
$MATLAB -nosplash -nojvm -nodisplay -nodesktop -r "try addpath('${util_path}'); \
    CBIG_pFIC_check_example_results('${user_output_path}', '${reference_output_path}'); catch; end; quit"

# remove output folder
rm -rf $CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/examples/output

exit 0
