#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
export MKL_NUM_THREADS=1
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run unit test
python -m unit_tests.CBIG_gcVAE_unit_test

# clean up
rm -rf $ROOTDIR"/unit_tests/checkpoints"
rm -rf $ROOTDIR"/unit_tests/data"
rm -rf $ROOTDIR"/unit_tests/results"
