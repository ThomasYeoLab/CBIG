#!/bin/sh

# unit test
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
export MKL_NUM_THREADS=1
mkdir -p $HOME/.cache/Python
export PYTHONPYCACHEPREFIX="${HOME}/.cache/Python"
source activate CBIG_An2024
module load cuda/11.0
cd $ROOTDIR

# run unit test
python -m unit_tests.test_CBIG_DeepResBat_unit_test

# clean up
rm -rf $ROOTDIR"/unit_tests/checkpoints"
rm -rf $ROOTDIR"/unit_tests/data"
rm -rf $ROOTDIR"/unit_tests/results"
