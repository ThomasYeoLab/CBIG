#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# setup
ROOTDIR=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/An2022_gcVAE"
log_dir=$ROOTDIR'/job_logs'
source activate CBIG_An2022
module load cuda/11.0
cd $ROOTDIR

# run
cmd_prefix="bash "$ROOTDIR"/replication/scripts/CBIG_gcVAE_replica_step6_train_cVAE.sh "$1
log_prefix="CBIG_gcVAE_replica_step6-1_"$2"perc_seed"
# seed0 - seed3
cmd_suffix0=" sample_size/"$2"perc_seed10"
cmd_suffix1=" sample_size/"$2"perc_seed11"
cmd_suffix2=" sample_size/"$2"perc_seed12"
cmd_suffix3=" sample_size/"$2"perc_seed13"
cmds_1=${cmd_prefix}${cmd_suffix0}";"
cmds=${cmds_1}${cmd_prefix}${cmd_suffix1}";"${cmd_prefix}${cmd_suffix2}";"${cmd_prefix}${cmd_suffix3}
log_suffix0="10"
log_suffix1="11"
log_suffix2="12"
log_suffix3="13"
logs_1=${log_prefix}${log_suffix0}";"
logs=${logs_1}${log_prefix}${log_suffix1}";"${log_prefix}${log_suffix2}";"${log_prefix}${log_suffix3}
python $ROOTDIR"/replication/scripts/"CBIG_gpuQ_parallel.py --cmds "$cmds" --log_names "$logs" --log_dir $log_dir
# seed4 - seed 7
cmd_suffix0=" sample_size/"$2"perc_seed14"
cmd_suffix1=" sample_size/"$2"perc_seed15"
cmd_suffix2=" sample_size/"$2"perc_seed16"
cmd_suffix3=" sample_size/"$2"perc_seed17"
cmds_1=${cmd_prefix}${cmd_suffix0}";"
cmds=${cmds_1}${cmd_prefix}${cmd_suffix1}";"${cmd_prefix}${cmd_suffix2}";"${cmd_prefix}${cmd_suffix3}
log_suffix0="14"
log_suffix1="15"
log_suffix2="16"
log_suffix3="17"
logs_1=${log_prefix}${log_suffix0}";"
logs=${logs_1}${log_prefix}${log_suffix1}";"${log_prefix}${log_suffix2}";"${log_prefix}${log_suffix3}
python $ROOTDIR"/replication/scripts/"CBIG_gpuQ_parallel.py --cmds "$cmds" --log_names "$logs" --log_dir $log_dir
# seed8 & seed 9
cmd_suffix0=" sample_size/"$2"perc_seed18"
cmd_suffix1=" sample_size/"$2"perc_seed19"
cmds=${cmd_prefix}${cmd_suffix0}";"${cmd_prefix}${cmd_suffix1}
log_suffix0="18"
log_suffix1="19"
logs=${log_prefix}${log_suffix0}";"${log_prefix}${log_suffix1}
python $ROOTDIR"/replication/scripts/"CBIG_gpuQ_parallel.py --cmds "$cmds" --log_names "$logs" --log_dir $log_dir
