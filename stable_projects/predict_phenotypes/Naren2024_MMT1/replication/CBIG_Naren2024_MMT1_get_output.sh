#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
UKBB_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1/UKBiobank"
RES_DIR="${BASE_DIR}/replication"
dir_name="WithinUKBB"
err_log="${LOG_DIR}/getoutput1_${dir_name}.err"
out_log="${LOG_DIR}/getoutput1_${dir_name}.out"
job_name="getoutput1_${dir_name}"

path_prefix1="/user_data/nwulan/ukbb_T1"
path_prefix2="/ukbb_train_test/files/from_hetong"
phe_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_test_final_9999_34.csv
sub_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_test_subj_list.txt
icv_dir=${UKBB_DIR}${path_prefix1}/data_prepare/ukbb_train_test/files/ICV/ukbb_all_sub_icv.csv
tra_sub_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_subj_list.txt
data_dir=${UKBB_DIR}${path_prefix1}/data_transfer/linear_output
out_dir="${RES_DIR}/results/${dir_name}"
out_subdir="basemodel_output"
model_dir="${RES_DIR}/results/${dir_name}/trained_model_ukbb"
phe_start_idx=1
phe_end_idx=35
batch_size=8
outputdim=33
index=98

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
icv_dir=$(readlink -f $icv_dir)
tra_sub_dir=$(readlink -f $tra_sub_dir)
data_dir=$(readlink -f $data_dir)

# get output used for meta-matching stacking for within UKBiobank experiment
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_test "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                            `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--outputdim ${outputdim}"



$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:20:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}


sleep 1s
# get output used for meta-matching finetune and classical transfer for within UKBiobank experiment
err_log="${LOG_DIR}/getoutput2_${dir_name}.err"
out_log="${LOG_DIR}/getoutput2_${dir_name}.out"
job_name="getoutput2_${dir_name}"
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_output "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                            `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx}"


$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:20:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}



#sleep 1s
# get output for HCPYA experiment
FILE_BASE_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
dir_name="HCPYA"
dataset=$dir_name
err_log="${LOG_DIR}/getoutput1_${dir_name}.err"
out_log="${LOG_DIR}/getoutput1_${dir_name}.out"
job_name="getoutput1_${dir_name}"

phe_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_final.csv
sub_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_subj_list.txt
icv_dir=${FILE_BASE_DIR}/${dataset}/ICV/icv.csv
ukbb_icv_dir=${UKBB_DIR}${path_prefix1}/data_prepare/ukbb_train_test/files/ICV/ukbb_all_sub_icv.csv
tra_sub_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_test_subj_list.txt
data_dir=${FILE_BASE_DIR}/${dataset}/T1_linear
out_dir="${RES_DIR}/results/${dir_name}"
model_dir="${RES_DIR}/results/UKBiobank67/trained_model_ukbb"
phe_start_idx=1
phe_end_idx=36
outputdim=67
index=98

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
ukbb_icv_dir=$(readlink -f $ukbb_icv_dir)
data_dir=$(readlink -f $data_dir)

# get output used for meta-matching stacking for HCPYA experiment
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_test "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                            `"--ukbb_icv_dir ${ukbb_icv_dir} "`
                                                           `"--outputdim ${outputdim} "`
                                                           `"--dataset ${dataset} "`
                                                           `"--across-dataset"


$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:05:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}


sleep 1s
# get output used for meta-matching finetune and classical transfer for HCPYA experiment
err_log="${LOG_DIR}/getoutput2_${dir_name}.err"
out_log="${LOG_DIR}/getoutput2_${dir_name}.out"
job_name="getoutput2_${dir_name}"
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_output "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                            `"--ukbb_icv_dir ${ukbb_icv_dir} "`
                                                           `"--dataset ${dataset} "`
                                                           `"--across-dataset"


$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:05:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}




sleep 1s
# get output for HCPA experiment
FILE_BASE_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
dir_name="HCPA"
dataset=$dir_name
err_log="${LOG_DIR}/getoutput1_${dir_name}.err"
out_log="${LOG_DIR}/getoutput1_${dir_name}.out"
job_name="getoutput1_${dir_name}"

phe_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_final.csv
sub_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_subj_list.txt
icv_dir=${FILE_BASE_DIR}/${dataset}/ICV/icv.csv
data_dir=${FILE_BASE_DIR}/${dataset}/T1_linear
out_dir="${RES_DIR}/results/${dir_name}"

phe_start_idx=1
phe_end_idx=46

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
data_dir=$(readlink -f $data_dir)

# get output used for meta-matching stacking for HCPA experiment
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_test "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                            `"--ukbb_icv_dir ${ukbb_icv_dir} "`
                                                           `"--outputdim ${outputdim} "`
                                                           `"--dataset ${dataset} "`
                                                           `"--across-dataset"


$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:05:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}


sleep 1s
# get output used for meta-matching finetune and classical transfer for HCPA experiment
err_log="${LOG_DIR}/getoutput2_${dir_name}.err"
out_log="${LOG_DIR}/getoutput2_${dir_name}.out"
job_name="getoutput2_${dir_name}"
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_dnn_output "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--ukbb_icv_dir ${ukbb_icv_dir} "`
                                                           `"--dataset ${dataset} "`
                                                           `"--across-dataset"


$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:05:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}
