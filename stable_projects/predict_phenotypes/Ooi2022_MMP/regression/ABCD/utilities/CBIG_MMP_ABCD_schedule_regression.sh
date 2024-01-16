#!/bin/sh
#####
# This function schedules the regressions in the CBIG scheduler. This script assumes that the first level, 
# second level and stats scripts are kept in the directory immediately above the one that this script is 
# kept in. User is to specify the type of regression and output directory.
# EXAMPLE: 
#    CBIG_MMP_ABCD_schedule_regression.sh KRR $parentdir/replication
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
if [ -z $CBIG_ABCD_REPDATA_DIR ]; then
    echo "ABCD_REPDATA_DIR not found: using default location on CBIG HPC"
    export CBIG_ABCD_REPDATA_DIR=$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Ooi2022_MMP/MMP_ABCD_data
fi
echo "ABCD_REPDATA_DIR = $CBIG_ABCD_REPDATA_DIR"
feature_folder=$CBIG_ABCD_REPDATA_DIR/data/features

### set up variables for scheduler based on type of regression
## first-level models
if [ $regression == "KRR" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_ABCD_KRR.sh "
    time_allowed="08:00:00"
    scheduler_tag="ABCD_KRR"
    mem="16GB"
elif [ $regression == "LRR" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_ABCD_LRR.sh "
    time_allowed="30:00:00"
    scheduler_tag="ABCD_LRR"
    mem="10GB"
elif [ $regression == "elasticnet" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_ABCD_Elasticnet.sh "
    time_allowed="30:00:00"
    scheduler_tag="ABCD_elasticnet"
    mem="10GB"

## second-level models
elif [ $regression == "multiKRR" ]; then
    base_cmd="$scriptdir/second_level/CBIG_MMP_ABCD_multiKRR.sh "
    time_allowed="48:00:00"
    scheduler_tag="ABCD_multiKRR"
    mem="20GB"

    # additional parameters needed for multiKRR
    # 4 kernel - all fmri
    outstem='multiKRR_k1_rs_k2_mid_k3_sst_k4_nback'
    feature_mat_cells="features_rs,features_mid,features_sst,features_nback"
    kernel_groups="1,2,3,4"
    mkrr_feature_folder=$feature_folder/fmri
elif [ $regression == "stacking" ]; then
    base_cmd="$scriptdir/second_level/CBIG_MMP_ABCD_stacking.sh "
    time_allowed="6:00:00"
    scheduler_tag="ABCD_stacking"
    mem="10GB"

    # additional parameters needed for stacking model
    second_lvl="LRR"
    outerFolds=120
    krr_output_dir=$outdir
    # models to be run
    outstem_arr=("stacking_${second_lvl}_rs_mid_sst_nback" \
            "stacking_${second_lvl}_best_cog" \
            "stacking_${second_lvl}_best_pers" \
            "stacking_${second_lvl}_best_mental" \
            "stacking_${second_lvl}_all")
    feature_mat_arr=("features_rs,features_mid,features_sst,features_nback" \
            "features_ct,features_tbss_ICVF,features_schaefer_FA,features_nback" \
            "features_ct,features_tbss_AD,features_schaefer_streamlen,features_nback" \
            "features_cv,features_tbss_ISOVF,features_schaefer_RD,features_sst" \
            "features_ct,features_cv,features_ca,features_rs,features_mid, \ 
            features_sst,features_nback,features_tbss_FA,features_tbss_MD, \
            features_tbss_AD,features_tbss_RD,features_tbss_OD,features_tbss_ICVF, \
            features_tbss_ISOVF,features_schaefer_streamcount_log,features_schaefer_streamlen, \
            features_schaefer_FA,features_schaefer_MD,features_schaefer_AD,features_schaefer_RD, \
            features_schaefer_OD,features_schaefer_ISOVF,features_schaefer_ICVF")

## stats
elif [ $regression == "stats" ]; then
    base_cmd="$scriptdir/stats/CBIG_MMP_ABCD_stats_perm_wrapper.sh "
    time_allowed="30:00:00"
    scheduler_tag="ABCD_stats"
    mem="10GB"

    # additional parameters needed for stacking model
    stats_outdir=$projectdir/stats
    run_once=0
    # multikrr outstem
    multiKRR_outstem='multiKRR_k1_rs_k2_mid_k3_sst_k4_nback'
    # stacking outstem and features
    outstem_arr=("stacking_LRR_rs_mid_sst_nback" "stacking_LRR_all")
    feature_mat_arr=("features_rs,features_mid,features_sst,features_nback" \
            "features_ct,features_cv,features_ca,features_rs,features_mid, \
            features_sst,features_nback,features_tbss_FA,features_tbss_MD, \
            features_tbss_AD,features_tbss_RD,features_tbss_OD,features_tbss_ICVF, \
            features_tbss_ISOVF,features_schaefer_streamcount_log,features_schaefer_streamlen, \
            features_schaefer_FA,features_schaefer_MD,features_schaefer_AD,features_schaefer_RD, \
            features_schaefer_OD,features_schaefer_ISOVF,features_schaefer_ICVF")
else
    echo "ERROR: regression option not recognised, exiting..."
    exit
fi

### set up feature names for regression and common files
# features
feature_array=("t1/features_ct" "t1/features_cv" "t1/features_ca" \
    "tbss/features_tbss_FA" "tbss/features_tbss_MD" "tbss/features_tbss_AD" \
    "tbss/features_tbss_RD" "tbss/features_tbss_OD" "tbss/features_tbss_ICVF" \
    "tbss/features_tbss_ISOVF" "tractography/features_schaefer_streamcount_log" \
    "tractography/features_schaefer_streamlen" "tractography/features_schaefer_FA" \
    "tractography/features_schaefer_MD" "tractography/features_schaefer_AD" \
    "tractography/features_schaefer_RD" "tractography/features_schaefer_OD" \
    "tractography/features_schaefer_ISOVF" "tractography/features_schaefer_ICVF" \
    "fmri/features_rs" "fmri/features_mid" "fmri/features_sst" "fmri/features_nback")

# common files and settings
sites=3
innerFolds=10
subtxt=$CBIG_ABCD_REPDATA_DIR/data/behaviors/MMP_ABCD_1823_subs.txt
subcsv=$CBIG_ABCD_REPDATA_DIR/data/behaviors/MMP_ABCD_behaviours.csv
predvar=$CBIG_ABCD_REPDATA_DIR/data/behaviors/MMP_ABCD_y_variables.txt
covtxt=$CBIG_ABCD_REPDATA_DIR/data/behaviors/MMP_ABCD_covariates.txt
ymat="MMP_ABCD_1823_y.mat"
covmat="MMP_ABCD_1823_covariates.mat"
y_output=$outdir/$ymat
cov_output=$outdir/$covmat

### define run job function: MODIFY THIS IF RUNNING ON YOUR OWN SCHEDULER
run_jobs(){
#run command
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime "$time_allowed" -name "$scheduler_tag" \
            -mem "$mem" -joberr "$logdir" -jobout "$logdir"

# wait for the files shared across random splits to be saved
if [ ! -f $y_output ] || [ ! -f $cov_output ]; then
    sleep 3m   
else
    sleep 2s
fi
}

### submit jobs to scheduler
for feature in ${feature_array[@]}; do
    if [ $regression == "KRR" ]; then
        feature_path=$feature_folder/$feature
        cmd="$base_cmd $feature_path $outdir $sites $innerFolds $subtxt $subcsv \
            $predvar $covtxt $ymat $covmat"
        run_jobs
    elif [ $regression == "LRR" ] || [ $regression == "elasticnet" ]; then
        for y_idx in $(seq 1 1 39); do 
            feature_path=$feature_folder/$feature
            cmd="$base_cmd $feature_path $outdir $sites $innerFolds $subtxt $subcsv \
                $predvar $covtxt $ymat $covmat $y_idx "
            run_jobs
        done
    elif [ $regression == "multiKRR" ]; then
        for fold_idx in $(seq 1 1 120); do 
            cmd="$base_cmd $outstem $mkrr_feature_folder $outdir $sites $innerFolds \
                $subtxt $subcsv $predvar $covtxt $ymat $covmat $fold_idx \
                $feature_mat_cells $kernel_groups "
            run_jobs
            # wait for kernels to be generated
            kernel_txt_file=$outdir/$outstem/results/step2_completion.txt
            while [ ! -f $kernel_txt_file ]; do
                echo "Waiting for kernels to be generated..."
                sleep 5m
            done
        done
        exit
    elif [ $regression == "stacking" ]; then
        arraylength=${#outstem_arr[@]}
        for (( i=0; i<${arraylength}; i++ )); do
            outstem=${outstem_arr[$i]}
            feature_mat_cells=$( echo ${feature_mat_arr[$i]} | xargs )
            for y_idx in $(seq 1 1 39); do 
                cmd="$base_cmd $outstem $krr_output_dir $outdir $outerFolds $innerFolds \
                $ymat $covmat $second_lvl $y_idx $feature_mat_cells "
                run_jobs
            done
        done
        exit
    elif [ $regression == "stats" ]; then
        for y_idx in $(seq 37 1 39); do
            # null models for single feature models
            featurebase=$(basename $feature)
            outstem="KRR_$featurebase"
            singlekrr_dir=$outdir/$outstem
            cmd="$base_cmd singleKRR $featurebase $stats_outdir $singlekrr_dir $subcsv $y_idx "
            run_jobs

            # null models for stacking models
            if [[ $run_once == 0 ]]; then
                arraylength=${#outstem_arr[@]}
                for (( i=0; i<${arraylength}; i++ )); do
                    stacking_outstem=${outstem_arr[$i]}
                    feature_mat_cells=$( echo ${feature_mat_arr[$i]} | xargs )
                    stacking_dir=$outdir/$stacking_outstem
                    cmd="$base_cmd stacking $stacking_outstem $stats_outdir $stacking_dir $subcsv $y_idx \
                    $outdir $feature_mat_cells "
                    run_jobs
                done
            fi
        done
        # null models for multiKRR models
        if [[ $run_once == 0 ]]; then
            multikrr_dir=$outdir/$multiKRR_outstem
            cmd="$base_cmd multiKRR $multiKRR_feature $stats_outdir $multikrr_dir $subcsv $y_idx "
            run_jobs
        fi
        run_once=1
    fi
done
