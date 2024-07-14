#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
err_log="${LOG_DIR}/unittest.err"
out_log="${LOG_DIR}/unittest.out"
job_name="MMT1_unittest"

cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m unit_tests.test_CBIG_MMT1_unit_test"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                   -walltime "00:15:00" \
                                   -mem "50G" -ncpus "4" -ngpus "1" -name $job_name \
                                   -joberr ${err_log} -jobout ${out_log}
