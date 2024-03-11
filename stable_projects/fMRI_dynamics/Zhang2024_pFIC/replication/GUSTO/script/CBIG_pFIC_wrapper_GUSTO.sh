#!/bin/bash
# This script can be used to replicate the pFIC model training, validaiton and test results based on the GUSTO
# dataset. This script does not require any input argument.
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# high-performance group
config=../config/GUSTO_high_performance.ini
module load cuda/11.7
source activate pFIC

# training
python ../../../model/CBIG_pFIC_training.py $config
python ../../../model/CBIG_pFIC_transpose_matrix.py $config

# validation
python ../../../model/CBIG_pFIC_validation.py $config

# testing
python ../../../model/CBIG_pFIC_test.py $config

source deactivate


# low-performance group
config=../config/GUSTO_low_performance.ini
source activate MFM

# training
python ../../../model/CBIG_pFIC_training.py $config
python ../../../model/CBIG_pFIC_transpose_matrix.py $config

# validation
python ../../../model/CBIG_pFIC_validation.py $config

# testing
python ../../../model/CBIG_pFIC_test.py $config

source deactivate


exit 0
