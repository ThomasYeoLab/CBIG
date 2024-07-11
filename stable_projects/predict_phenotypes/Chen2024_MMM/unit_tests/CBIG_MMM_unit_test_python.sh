#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

ut_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/unit_tests/"
log_dir="${ut_dir}/log"
mkdir -p ${log_dir}
cmd="cd ${ut_dir}; conda activate CBIG_Chen2024;"
cmd="${cmd} export PYTHONPATH=${PYTHONPATH}:$PWD/../;"
cmd="${cmd} export MKL_SERVICE_FORCE_INTEL=1;"
cmd="${cmd} python test_CBIG_MMM_unit_test.py;"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 1:00:00 -mem 32G -ncpus 4 -ngpus 1 \
-name "unit_test_py" -joberr "${log_dir}/unit_test_py_err.txt" \
-jobout "${log_dir}/unit_test_py_out.txt"
