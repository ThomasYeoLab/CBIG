#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
BASE_DIR="${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
LOG_DIR="${BASE_DIR}/job_logs"
UKBB_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1/UKBiobank"
FILE_BASE_DIR="${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/Naren2024_MMT1"
RES_DIR="${BASE_DIR}/replication"
dir_name="HCPA"
dataset=$dir_name
err_log="${LOG_DIR}/elasticnet_${dir_name}.err"
out_log="${LOG_DIR}/elasticnet_${dir_name}.out"
job_name="elasticnet_${dir_name}"

phe_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_final.csv
sub_dir=${FILE_BASE_DIR}/${dataset}/${dataset}_diff_roi_subj_list.txt
data_dir=${FILE_BASE_DIR}/${dataset}/IDPs/${dataset}_idp_ASEG_DKT_update.csv
out_dir="${RES_DIR}/results/${dir_name}"
inter_dir="${out_dir}/output_intermediate"
ks=(10 20 50 100 200)
rng=100
phe_start_idx=1
phe_end_idx=46

# change symbolic link to real path
phe_dir=$(readlink -f $phe_dir)
sub_dir=$(readlink -f $sub_dir)
data_dir=$(readlink -f $data_dir)

## train for elasticnet
err_log="${LOG_DIR}/elasticnet_${dir_name}.err"
out_log="${LOG_DIR}/elasticnet_${dir_name}.out"
job_name="elasticnet_${dir_name}"
out_subdir="/elasticnet"
# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m ElasticNet.elasticnet "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--rng ${rng} "`
                                                           `"--ks ${ks[@]} --across-dataset "`
                                                           `"--dataset ${dir_name}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "8:00:00" \
                                     -mem "128G" -ncpus "4" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}


sleep 1s
# train for classical transfer learning
err_log="${LOG_DIR}/classical_trans_${dir_name}.err"
out_log="${LOG_DIR}/classical_trans_${dir_name}.out"
job_name="classical_trans_${dir_name}"

path_prefix1="/user_data/nwulan/ukbb_T1"
path_prefix2="/ukbb_train_test/files/from_hetong"
icv_dir=${FILE_BASE_DIR}/${dataset}/ICV/icv.csv
ukbb_icv_dir=${UKBB_DIR}${path_prefix1}/data_prepare/ukbb_train_test/files/ICV/ukbb_all_sub_icv.csv
tra_sub_dir=${UKBB_DIR}${path_prefix1}/data_prepare${path_prefix2}/ukbb_train_test_subj_list.txt
data_dir="${RES_DIR}/results/${dir_name}/basemodel_output"
out_dir="${RES_DIR}/results/${dir_name}"
out_subdir="/classical_transfer"
model_dir="${RES_DIR}/results/UKBiobank67/trained_model_ukbb"
in_c=256
out_c=64
epochs=50
batch_size=200
index=98

# change symbolic link to real path
icv_dir=$(readlink -f $icv_dir)
ukbb_icv_dir=$(readlink -f $ukbb_icv_dir)
tra_sub_dir=$(readlink -f $tra_sub_dir)

# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_classical_transfer "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--epochs ${epochs} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--ks ${ks[@]} --across-dataset "`
                                                           `"--n_rng ${rng} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                           `"--ukbb_icv_dir ${ukbb_icv_dir} "`
                                                             `"--in_c ${in_c} "`
                                                           `"--out_c ${out_c} "`
                                                           `"--dataset ${dataset}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "120:00:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}



sleep 1s
# train for meta-matching finetune
err_log="${LOG_DIR}/mm_finetune_${dir_name}.err"
out_log="${LOG_DIR}/mm_finetune_${dir_name}.out"
job_name="mm_finetune_${dir_name}"
out_subdir="/mm_finetune"
metric='cod'

# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_mm_finetune "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--epochs ${epochs} "`
                                                           `"--batch_size ${batch_size} "`
                                                           `"--ks ${ks[@]} --across-dataset "`
                                                           `"--n_rng ${rng} "`
                                                           `"--model_dir ${model_dir} "`
                                                           `"--index ${index} "`
                                                           `"--metric ${metric} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--icv_dir ${icv_dir} "`
                                                           `"--tra_sub_dir ${tra_sub_dir} "`
                                                           `"--ukbb_icv_dir ${ukbb_icv_dir} "`
                                                             `"--in_c ${in_c} "`
                                                           `"--out_c ${out_c} "`
                                                           `"--dataset ${dataset}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "120:00:00" \
                                     -mem "128G" -ncpus "4" -ngpus "1" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}



sleep 1s
# train for meta-matching stacking
err_log="${LOG_DIR}/mm_stacking_${dir_name}.err"
out_log="${LOG_DIR}/mm_stacking_${dir_name}.out"
job_name="mm_stacking_${dir_name}"
out_subdir="/mm_stacking"
k_limit=33

# shellcheck disable=SC2124
cmd="source activate Naren2024_MMT1;cd ${BASE_DIR};python -m cbig.CBIG_ukbb_mm_stacking "`
                                                           `"--phe_dir ${phe_dir} "`
                                                           `"--sub_dir ${sub_dir} "`
                                                           `"--data_dir ${data_dir} "`
                                                           `"--out_dir ${out_dir} "`
                                                           `"--out_subdir ${out_subdir} "`
                                                           `"--inter_dir ${inter_dir} "`
                                                           `"--ks ${ks[@]} --across-dataset "`
                                                           `"--n_rng ${rng} "`
                                                           `"--metric ${metric} "`
                                                           `"--start_idx ${phe_start_idx} "`
                                                           `"--end_idx ${phe_end_idx} "`
                                                           `"--k_limit ${k_limit} "`
                                                           `"--dataset ${dataset}"

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
                                     -walltime "00:30:00" \
                                     -mem "128G" -ncpus "4" -name $job_name \
                                     -joberr ${err_log} -jobout ${out_log}