#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
exp_names=("WithinUKBB" "HCPYA" "HCPA")

for exp_name_id in {0..2}; do
    exp_name=${exp_names[exp_name_id]}
    err_log="${LOG_DIR}/check_replication_res_${exp_name}.err"
    out_log="${LOG_DIR}/check_replication_res_${exp_name}.out"
    job_name="MMT1_CheckRepliacation_${exp_name}"

    cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m utils.check_reference_results "`
                                                                   `"--datasetname ${exp_name} "

    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                       -walltime "00:10:00" \
                                       -mem "10G" -ncpus "1" -name $job_name \
                                       -joberr ${err_log} -jobout ${out_log}
    sleep 1s
done
