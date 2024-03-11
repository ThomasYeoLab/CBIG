#!/bin/bash

# This wrapper script can be used to replicate the analysis results from HCP dataset.
# This script does not require any input argument.
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

MATLAB=`which matlab`

config=../config/HCP.ini
module load cuda/11.7
source activate pFIC

# training
python -W ignore ../../../model/CBIG_pFIC_training.py $config
python -W ignore ../../../model/CBIG_pFIC_transpose_matrix.py $config

# validation
python -W ignore ../../../model/CBIG_pFIC_validation.py $config

# testing
python -W ignore ../../../model/CBIG_pFIC_test.py $config

source deactivate

# check results against reference output
util_path=`realpath ../script`
user_output_path=`realpath ../test`
reference_output_path=`realpath ../reference_output/test`
$MATLAB -nosplash -nojvm -nodisplay -nodesktop -r "try addpath('${util_path}'); \
    CBIG_pFIC_check_HCP_results('${user_output_path}', '${reference_output_path}'); catch; end; quit"

exit 0
