#!/bin/bash

# submission
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
log_dir=$ROOTDIR'/job_logs'
job_prefix="DRB"
cd $ROOTDIR

# submission
cmd_prefix="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step6_eval_dataset_pred.sh "$1

# unharm
cmd=$cmd_prefix" unharm"
job_name=$job_prefix"_Step6_"$1"_unharm"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
    -name $job_name \
    -walltime 04:00:00 \
    -mem 4G \
    -ncpus 1 \
    -jobout $log_dir"/"$job_name".out" \
    -joberr $log_dir"/"$job_name".err"
sleep 10s
# ComBat
cmd=$cmd_prefix" ComBat4cov"
job_name=$job_prefix"_Step6_"$1"_ComBat4cov"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
    -name $job_name \
    -walltime 04:00:00 \
    -mem 4G \
    -ncpus 1 \
    -jobout $log_dir"/"$job_name".out" \
    -joberr $log_dir"/"$job_name".err"
sleep 10s
# CovBat
cmd=$cmd_prefix" CovBat4cov"
job_name=$job_prefix"_Step6_"$1"_CovBat4cov"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
    -name $job_name \
    -walltime 04:00:00 \
    -mem 4G \
    -ncpus 1 \
    -jobout $log_dir"/"$job_name".out" \
    -joberr $log_dir"/"$job_name".err"
sleep 10s
# cVAE
cmd=$cmd_prefix" cVAE"
job_name=$job_prefix"_Step6_"$1"_cVAE"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
    -name $job_name \
    -walltime 04:00:00 \
    -mem 4G \
    -ncpus 1 \
    -jobout $log_dir"/"$job_name".out" \
    -joberr $log_dir"/"$job_name".err"
sleep 10s
# coVAE
cmd=$cmd_prefix" coVAE"
job_name=$job_prefix"_Step6_"$1"_coVAE"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
    -name $job_name \
    -walltime 04:00:00 \
    -mem 4G \
    -ncpus 1 \
    -jobout $log_dir"/"$job_name".out" \
    -joberr $log_dir"/"$job_name".err"
sleep 10s
# DeepResBat
cmd=$cmd_prefix" DeepResBat"
job_name=$job_prefix"_Step6_"$1"_DeepResBat"
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "$cmd" \
    -name $job_name \
    -walltime 04:00:00 \
    -mem 4G \
    -ncpus 1 \
    -jobout $log_dir"/"$job_name".out" \
    -joberr $log_dir"/"$job_name".err"
