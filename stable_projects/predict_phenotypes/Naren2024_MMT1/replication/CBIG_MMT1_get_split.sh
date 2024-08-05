#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
FILE_BASE_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
RES_DIR="${BASE_DIR}/replication"
dir_name="WithinUKBB"
dataset="ukbb"
err_log="${LOG_DIR}/get_split_${dir_name}.err"
out_log="${LOG_DIR}/get_split_${dir_name}.out"
job_name="get_split_${dir_name}"

path_prefix1="/user_data/nwulan/ukbb_T1"
path_prefix2="/ukbb_train_test/files/from_hetong"
phe_dir=${FILE_BASE_DIR}/UKBiobank/${path_prefix1}/data_prepare${path_prefix2}/ukbb_test_final_9999_34.csv
sub_dir=${FILE_BASE_DIR}/UKBiobank/${path_prefix1}/data_prepare${path_prefix2}/ukbb_test_subj_list.txt
out_dir="${RES_DIR}/results/${dir_name}"
inter_dir="${out_dir}/output_intermediate"

ks=(10 20 50 100 200)
rng=100
num_inner_folds=5
seed=0
phe_start_idx=1
phe_end_idx=35
n_sub_s=9888


mkdir -p $LOG_DIR

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)

# get training and test split for k-shot for within UKBiobank experiment
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m utils.generate_split_kshot "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--n_rng ${rng} "`
                                                           `"--ks ${ks[@]} "`
                                                           `"--dataset ${dataset} "`
                                                           `"--num_inner_folds ${num_inner_folds} "`
                                                           `"--seed ${seed}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:10:00" \
                                     -mem "50G" -ncpus "4" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}


#sleep 1s
## get training and test split for k-shot for HCPYA experiment
#dir_name="HCPYA"
#dataset=$dir_name
#err_log="${LOG_DIR}/get_split_${dir_name}.err"
#out_log="${LOG_DIR}/get_split_${dir_name}.out"
#job_name="get_split_${dir_name}"
#phe_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_final.csv
#sub_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_subj_list.txt
#out_dir="${RES_DIR}/results/${dir_name}"
#inter_dir="${out_dir}/output_intermediate"
#
#phe_start_idx=1
#phe_end_idx=36
#n_sub_s=1017
#
## change symbolic link to real path
#phe_dir=$(readlink -f $phe_dir)
#sub_dir=$(readlink -f $sub_dir)
#
## shellcheck disable=SC2124
#cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m utils.generate_split_kshot "`
#                                                           `"--phe_dir ${phe_dir} "`
#                                                           `"--sub_dir ${sub_dir} "`
#                                                           `"--inter_dir ${inter_dir} "`
#                                                           `"--start_idx ${phe_start_idx} "`
#                                                           `"--end_idx ${phe_end_idx} "`
#                                                           `"--n_rng ${rng} "`
#                                                           `"--ks ${ks[@]} "`
#                                                           `"--dataset ${dir_name} "`
#                                                           `"--num_inner_folds ${num_inner_folds} "`
#                                                           `"--seed ${seed}"
#
#$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
#                                     -walltime "00:10:00" \
#                                     -mem "50G" -ncpus "4" -name $job_name \
#                                     -joberr ${err_log} -jobout ${out_log}
#
#
#sleep 1s
## get training and test split for k-shot for HCPA experiment
#dir_name="HCPA"
#dataset=$dir_name
#err_log="${LOG_DIR}/get_split_${dir_name}.err"
#out_log="${LOG_DIR}/get_split_${dir_name}.out"
#job_name="get_split_${dir_name}"
#phe_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_final.csv
#sub_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_subj_list.txt
#out_dir="${RES_DIR}/results/${dir_name}"
#inter_dir="${out_dir}/output_intermediate"
#
#phe_start_idx=1
#phe_end_idx=46
#n_sub_s=656
#
## change symbolic link to real path
#phe_dir=$(readlink -f $phe_dir)
#sub_dir=$(readlink -f $sub_dir)
#
## shellcheck disable=SC2124
#cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m utils.generate_split_kshot "`
#                                                           `"--phe_dir ${phe_dir} "`
#                                                           `"--sub_dir ${sub_dir} "`
#                                                           `"--inter_dir ${inter_dir} "`
#                                                           `"--start_idx ${phe_start_idx} "`
#                                                           `"--end_idx ${phe_end_idx} "`
#                                                           `"--n_rng ${rng} "`
#                                                           `"--ks ${ks[@]} "`
#                                                           `"--dataset ${dir_name} "`
#                                                           `"--num_inner_folds ${num_inner_folds} "`
#                                                           `"--seed ${seed}"
#
#$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
#                                     -walltime "00:10:00" \
#                                     -mem "50G" -ncpus "4" -name $job_name \
#                                     -joberr ${err_log} -jobout ${out_log}