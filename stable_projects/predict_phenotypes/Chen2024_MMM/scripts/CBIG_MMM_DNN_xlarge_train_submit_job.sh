#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

rep_dir="$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Chen2024_MMM/replication/"
log_dir="${rep_dir}/log"
mkdir -p ${log_dir}
cmd="cd ${rep_dir}; source activate CBIG_Chen2024;"
cmd="${cmd} python ../cbig/Chen2024/CBIG_dnn_xlarge_train.py --src_dataset $1 --gpu 0 --seed 1 --epochs 100 --metric cod \
--weight_decay 2.727460424379488e-07 --lr 0.036612177992895435 --dropout 0.4 \
--n_l1 512 --n_l2 256 --n_l3 128 --n_l4 1024 --n_hidden_layer 2  --batch_size 128 --patience 25"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" -walltime 3:00:00 -mem 192G -ncpus 2 -ngpus 1 -name "MMM_DNN_train" \
-joberr "${log_dir}/DNN_train_err.txt" -jobout "${log_dir}/DNN_train_out.txt"
