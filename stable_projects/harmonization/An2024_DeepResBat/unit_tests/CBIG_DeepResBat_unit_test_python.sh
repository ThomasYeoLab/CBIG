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
test_status=$?

# clean up generated output only if the test passed; the test's own tearDown
# already does this, this is just a safety net (e.g. in case of a hard crash
# that skips tearDown). Output is kept on failure for debugging.
if [ $test_status -eq 0 ]; then
    rm -rf $ROOTDIR"/unit_tests/output"
fi

exit $test_status
