#!/bin/sh
#####
# This function schedules the regressions in the CBIG scheduler. This script assumes that the first level, 
# second level and stats scripts are kept in the directory immediately above the one that this script is 
# kept in. User is to specify the type of regression and output directory.
# EXAMPLE: 
#    CBIG_MMP_HCP_schedule_regression_example.sh KRR $parentdir/replication
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

### choose type of regression and specify output directory
regression=$1
projectdir=$2
outdir=$projectdir/output
logdir=$projectdir/scheduler_logs
if [ ! -d $logdir ]; then mkdir -p $logdir; fi
scriptdir=$(dirname $(dirname "$(readlink -f "$0")"))
if [ -z $CBIG_HCP_REPDATA_DIR ]; then
    echo "HCP_REPDATA_DIR not found: using default location on CBIG HPC"
    export CBIG_HCP_REPDATA_DIR=$CBIG_CODE_DIR/stable_projects/predict_phenotypes/Ooi2022_MMP/examples/input
fi
echo "HCP_REPDATA_DIR = $CBIG_HCP_REPDATA_DIR"
feature_folder=$CBIG_HCP_REPDATA_DIR/features

### set up variables for scheduler based on type of regression
## first-level models
if [ $regression == "KRR" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_HCP_KRR_example.sh "
    time_allowed="08:00:00"
    scheduler_tag="Ooi2022_Example"
    mem="10GB"
elif [ $regression == "LRR" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_HCP_LRR_example.sh "
    time_allowed="60:00:00"
    scheduler_tag="Ooi2022_Example"
    mem="10GB"
elif [ $regression == "elasticnet" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_HCP_Elasticnet_example.sh "
    time_allowed="60:00:00"
    scheduler_tag="Ooi2022_Example"
    mem="10GB"

## second-level models
elif [ $regression == "multiKRR" ]; then
    base_cmd="$scriptdir/second_level/CBIG_MMP_HCP_multiKRR_example.sh "
    time_allowed="16:00:00"
    scheduler_tag="Ooi2022_Example"
    mem="20GB"

    # additional parameters needed for multiKRR
    # 6 kernel - all fmri
    outstem='multiKRR_k1_cv_k2_rs'
    feature_mat_cells="features_cv,features_rs"
    kernel_groups="1,2"
    mkrr_feature_folder=$feature_folder/fmri
elif [ $regression == "stacking" ]; then
    base_cmd="$scriptdir/second_level/CBIG_MMP_HCP_stacking_example.sh "
    time_allowed="6:00:00"
    scheduler_tag="Ooi2022_Example"
    mem="10GB"

    # additional parameters needed for stacking model
    second_lvl="LRR"
    outerFolds=120
    krr_output_dir=$outdir
    # models to be run
    outstem_arr=("stacking_${second_lvl}_cv_rs")
    feature_mat_arr=("features_cv,features_rs")
fi

# set up features to run
feature_array=("t1/features_cv" "fmri/features_rs") 

# common files and settings
outerFolds=3
innerFolds=3
subtxt=$CBIG_HCP_REPDATA_DIR/behaviors/MMP_HCP_60_subs.txt
scorecsv=$CBIG_HCP_REPDATA_DIR/behaviors/MMP_HCP_componentscores.csv
# please modify path to restricted CSV: NOT USED IN EXAMPLE
restricted_dir=$CBIG_REPDATA_DIR/some_directory
restrictedcsv=$restricted_dir/some_csv.csv
predvar=$CBIG_HCP_REPDATA_DIR/behaviors/MMP_HCP_example_y_variables.txt
covtxt=$CBIG_HCP_REPDATA_DIR/behaviors/MMP_HCP_example_covariates.txt
ymat="MMP_HCP_60_y.mat"
covmat="MMP_HCP_60_covariates.mat"
y_output=$outdir/$ymat
cov_output=$outdir/$covmat

### define run job function: MODIFY THIS IF RUNNING ON YOUR OWN SCHEDULER
run_jobs(){
#run command
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime "$time_allowed" -name "$scheduler_tag" \
            -mem "$mem" -joberr "$logdir" -jobout "$logdir"

# wait for the files shared across random splits to be saved
if [ ! -f $y_output ] || [ ! -f $cov_output ]; then
    sleep 1m   
else
    sleep 2s
fi
}

### submit jobs to scheduler
for feature in ${feature_array[@]}; do
    if [ $regression == "KRR" ]; then
        feature_path=$feature_folder/$feature
        for split_idx in $(seq 1 1 1); do 
            cmd="$base_cmd $feature_path $outdir $outerFolds $innerFolds $subtxt $scorecsv \
            $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx "
            run_jobs
        done
    elif [ $regression == "LRR" ] || [ $regression == "elasticnet" ]; then
        feature_path=$feature_folder/$feature
        num_var=3
        for split_idx in $(seq 1 1 1); do 
            cmd="$base_cmd $feature_path $outdir $outerFolds $innerFolds $subtxt $scorecsv \
            $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx $num_var "
            run_jobs
        done
    elif [ $regression == "multiKRR" ]; then
        # run jobs to generate kernels for all splits first
        for split_idx in $(seq 1 1 1); do 
            cmd="$base_cmd $outstem $mkrr_feature_folder $outdir $outerFolds $innerFolds \
            $subtxt $scorecsv $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx 1 \
            $feature_mat_cells $kernel_groups "
            run_jobs
        done
        # remaining jobs only submitted if kernel for that split is already generated
        for split_idx in $(seq 1 1 1); do 
            for fold_idx in $(seq 2 1 10); do 
                # wait for kernels to be generated
                kernel_txt_file=$outdir/$outstem/seed_${split_idx}/results/step2_completion.txt
                while [ ! -f $kernel_txt_file ]; do
                    echo "Waiting for kernels to be generated..."
                    sleep 1m
                done
                cmd="$base_cmd $outstem $mkrr_feature_folder $outdir $outerFolds $innerFolds \
                    $subtxt $scorecsv $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx \
                    $fold_idx $feature_mat_cells $kernel_groups "
                run_jobs
            done
        done
        exit
    elif [ $regression == "stacking" ]; then
        num_var=3
        arraylength=${#outstem_arr[@]}
        for (( i=0; i<${arraylength}; i++ )); do
            outstem=${outstem_arr[$i]}
            feature_mat_cells=$( echo ${feature_mat_arr[$i]} | xargs )
            for split_idx in $(seq 1 1 1); do
                cmd="$base_cmd $outstem $krr_output_dir $outdir $outerFolds $innerFolds \
                $ymat $covmat $second_lvl $split_idx $num_var $feature_mat_cells "
                run_jobs
            done
        done
        exit
    fi
done
