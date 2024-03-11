#!/bin/bash

# This wrapper script provides an example for parameteric Feedback Inhibition Control (pFIC) model training.
# This script also serves as a unit-test test case, so it automatically checks the generated output with 
# a set of reference output.

# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

MATLAB=`which matlab`

config=../config/example.ini
module load cuda/11.7
source activate pFIC

# training
python ../../model/CBIG_pFIC_training.py $config
python ../../model/CBIG_pFIC_transpose_matrix.py $config

conda deactivate

# check results against reference output
util_path=`realpath ../../examples`
user_output_path=`realpath ../../examples/output`
reference_output_path=`realpath ../../examples/reference_output`
$MATLAB -nosplash -nojvm -nodisplay -nodesktop -r "try addpath('${util_path}'); \
    CBIG_pFIC_check_example_results('${user_output_path}', '${reference_output_path}'); catch; end; quit"

exit 0
