#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
FILE_BASE_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
RES_DIR="${BASE_DIR}/replication"
dir_name="HCPA"
dataset=$dir_name

err_log="${LOG_DIR}/haufe_check_${dir_name}.err"
out_log="${LOG_DIR}/haufe_check_${dir_name}.out"
job_name="haufe_check_${dir_name}"

sub_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_subj_list.txt
data_dir=${FILE_BASE_DIR}/${dataset}/T1_linear
out_dir="${RES_DIR}/results/${dir_name}"
inter_dir="${out_dir}/output_intermediate"
ik=3
k=100
rng=100
stem="mm_stacking"
ib=44
p="moca_total"

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
data_dir=$(readlink -f $data_dir)

# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_haufe "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--ib ${ib} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--stem ${stem} "`
                                                           `"--phe ${p} "`
                                                           `"--rng ${rng} "`
                                                           `"--ik ${ik} "`
                                                           `"--k ${k} "`
                                                           `"--dataset ${dir_name}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "1:00:00" \
                                     -mem "128G" -ncpus "4" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}



