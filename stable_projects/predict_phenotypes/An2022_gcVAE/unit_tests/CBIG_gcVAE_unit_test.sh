#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
export MKL_NUM_THREADS=1
source CBIG_init_conda
conda activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run unit test
python -m unit_tests.test_CBIG_gcVAE_unit_test
test_status=$?

# clean up generated output only if the test passed; kept on failure for
# debugging.
if [ $test_status -eq 0 ]; then
    rm -rf $ROOTDIR"/unit_tests/output"
fi

exit $test_status
