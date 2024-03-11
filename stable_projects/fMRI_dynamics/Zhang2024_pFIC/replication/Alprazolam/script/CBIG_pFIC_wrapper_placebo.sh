#!/bin/bash

# This script can be used to replicate the pFIC model training, validaiton and test results based on the Alprazolam
# placebo session data. This script does not require any input argument.
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

MATLAB=`which matlab`

# placebo
config=../config/placebo.ini
module load cuda/11.7
source activate pFIC

# training
python ../../../model/CBIG_pFIC_training.py $config
python ../../../model/CBIG_pFIC_transpose_matrix.py $config

# validation
python ../../../model/CBIG_pFIC_validation.py $config

# testing
python ../../../model/CBIG_pFIC_test.py $config

# generate extrapolated model parameters
util_path=`realpath ../../../util/Alprazolam`
test_file_path=`realpath ../reference_output/test/placebo/test_all.csv`
training_folder_path=`realpath ../reference_output/training/placebo`
training_file_name=training
output_path=`realpath ../reference_output/test/placebo`
myelin_full=`realpath ../../HCP/input/myelin.csv`
rsfc_gradient_full=`realpath ../../HCP/input/rsfc_gradient.csv`
$MATLAB -nosplash -nojvm -nodisplay -nodesktop -r "try addpath('${util_path}'); \
    CBIG_pFIC_extrapolate_best_parameter('${test_file_path}', '${training_folder_path}', \
    '${training_file_name}', '${output_path}', '${myelin_full}', '${rsfc_gradient_full}'); \
    catch; end; quit"

# run extrapolation
python ../../../model/CBIG_pFIC_extrapolation.py $config

source deactivate



exit 0
