#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
UKBB_DIR=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1/UKBiobank
RES_DIR="${BASE_DIR}/replication"
dir_name="WithinUKBB"
err_log="${LOG_DIR}/training_${dir_name}.err"
out_log="${LOG_DIR}/training_${dir_name}.out"
job_name="training_${dir_name}"

path_prefix1="/user_data/nwulan/ukbb_T1"
path_prefix2="/ukbb_train_test/files/from_hetong"
phe_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_test_final_36838_33.csv
icv_dir=${UKBB_DIR}${path_prefix1}/data_prepare/ukbb_train_test/files/ICV/ukbb_all_sub_icv.csv
sub_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_subj_list.txt
data_dir=${UKBB_DIR}${path_prefix1}/data_transfer/linear_output
out_dir="${RES_DIR}/results/${dir_name}"
out_subdir="trained_model_ukbb"
inter_dir="${out_dir}/output_intermediate"
epochs=100
batch_size=8
lr=1e-5
scheduler_decrease=50
cn=("32 64 128 256 256 64")
output_dim=33
dropout=0.1
phe_start_idx=1
phe_end_idx=34


# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
icv_dir=$(readlink -f $icv_dir)
data_dir=$(readlink -f $data_dir)

# train for Within UKBB experiments
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_train "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--epochs ${epochs} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--lr ${lr} "`
                                                           `"--scheduler_decrease ${scheduler_decrease} "`
                                                           `"--channel_number ${cn[@]} "`
                                                           `"--output_dim ${output_dim} "`
                                                           `"--dropout ${dropout} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--icv_dir ${icv_dir}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "120:00:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}


sleep 1s
# train for across dataset experiments
dir_name="UKBiobank67"
err_log="${LOG_DIR}/training_${dir_name}.err"
out_log="${LOG_DIR}/training_${dir_name}.out"
job_name="training_${dir_name}"

phe_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_test_final.csv
icv_dir=${UKBB_DIR}${path_prefix1}/data_prepare/ukbb_train_test/files/ICV/ukbb_all_sub_icv.csv
sub_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_test_subj_list.txt
out_dir="${RES_DIR}/results/${dir_name}"
output_dim=67
phe_end_idx=68

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
icv_dir=$(readlink -f $icv_dir)

# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_train "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--epochs ${epochs} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--lr ${lr} "`
                                                           `"--scheduler_decrease ${scheduler_decrease} "`
                                                           `"--channel_number ${cn[@]} "`
                                                           `"--output_dim ${output_dim} "`
                                                           `"--dropout ${dropout} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--icv_dir ${icv_dir}"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "120:00:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}