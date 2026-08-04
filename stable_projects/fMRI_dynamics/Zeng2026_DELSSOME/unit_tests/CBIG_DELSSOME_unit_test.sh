#!/bin/sh
# Written by Tianchu Zeng and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Wrapper that runs the Zeng2026_DELSSOME group-level unit test. Assumes the
# CBIG_DELSSOME conda environment (CBIG_DELSSOME_python_env.yml, next to this file).
# The predictor-training and DELSSOME examples also need the larger training matrices
# / pretrained predictor; see test_CBIG_DELSSOME_unit_test.py and README.md for how
# those are resolved.

HERE=$(dirname "$(readlink -f "$0")")

source activate CBIG_DELSSOME
python -u "${HERE}/test_CBIG_DELSSOME_unit_test.py"
exit_code=$?
conda deactivate

exit ${exit_code}
