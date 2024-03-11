#!/bin/bash

# This wrapper script can be used to replicate the analysis results from HCP dataset where the parameters
# are only parameterized by RSFC gradient. This script does not require any input argument.
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

config=../config/HCP_rsfc_only.ini
module load cuda/11.7
source activate pFIC

# training
python -W ignore CBIG_pFIC_training_HCP_rsfc_gradient_only.py $config
python -W ignore CBIG_pFIC_transpose_matrix.py $config

# validation
python -W ignore CBIG_pFIC_validation_HCP_rsfc_gradient_only.py $config

# testing
python -W ignore CBIG_pFIC_test_HCP.py $config

source deactivate

# check results against reference output
util_path=`realpath ../script`
user_output_path=`realpath ../test_rsfc_gradient_only`
reference_output_path=`realpath ../reference_output/test_rsfc_gradient_only`
$MATLAB -nosplash -nojvm -nodisplay -nodesktop -r "try addpath('${util_path}'); \
    CBIG_pFIC_check_HCP_results('${user_output_path}', '${reference_output_path}'); catch; end; quit"

exit 0
