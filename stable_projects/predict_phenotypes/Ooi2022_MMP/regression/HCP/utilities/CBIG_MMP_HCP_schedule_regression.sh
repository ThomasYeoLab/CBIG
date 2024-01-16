#!/bin/sh
#####
# This function schedules the regressions in the CBIG scheduler. This script assumes that the first level, 
# second level and stats scripts are kept in the directory immediately above the one that this script is 
# kept in. User is to specify the type of regression and output directory.
# EXAMPLE: 
#    CBIG_MMP_HCP_schedule_regression.sh KRR $parentdir/replication
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
    export CBIG_HCP_REPDATA_DIR=$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Ooi2022_MMP/MMP_HCP_data
fi
echo "HCP_REPDATA_DIR = $CBIG_HCP_REPDATA_DIR"
behav_dir=$CBIG_HCP_REPDATA_DIR/data/behaviors
feature_folder=$CBIG_HCP_REPDATA_DIR/data/features

### set up variables for scheduler based on type of regression
## first-level models
if [ $regression == "KRR" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_HCP_KRR.sh "
    time_allowed="08:00:00"
    scheduler_tag="HCP_KRR"
    mem="10GB"
elif [ $regression == "LRR" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_HCP_LRR.sh "
    time_allowed="60:00:00"
    scheduler_tag="HCP_LRR"
    mem="10GB"
elif [ $regression == "elasticnet" ]; then
    base_cmd="$scriptdir/first_level/CBIG_MMP_HCP_Elasticnet.sh "
    time_allowed="60:00:00"
    scheduler_tag="HCP_elasticnet"
    mem="10GB"

## second-level models
elif [ $regression == "multiKRR" ]; then
    base_cmd="$scriptdir/second_level/CBIG_MMP_HCP_multiKRR.sh "
    time_allowed="16:00:00"
    scheduler_tag="HCP_multiKRR"
    mem="20GB"

    # additional parameters needed for multiKRR
    # 6 kernel - all fmri
    outstem='multiKRR_k1_rs_k2_wm_k3_lang_k4_gamb_k5_social_k6_motor'
    feature_mat_cells="features_rs,features_wm,features_lang,features_gamb,features_social,features_motor"
    kernel_groups="1,2,3,4,5,6"
    mkrr_feature_folder=$feature_folder/fmri
elif [ $regression == "stacking" ]; then
    base_cmd="$scriptdir/second_level/CBIG_MMP_HCP_stacking.sh "
    time_allowed="6:00:00"
    scheduler_tag="HCP_stacking"
    mem="10GB"

    # additional parameters needed for stacking model
    second_lvl="LRR"
    outerFolds=120
    krr_output_dir=$outdir
    # models to be run
    outstem_arr=("stacking_${second_lvl}_rs_wm_lang_gamb_social_motor" \
            "stacking_${second_lvl}_best_cog" \
            "stacking_${second_lvl}_best_satisf" \
            "stacking_${second_lvl}_best_er" \
            "stacking_${second_lvl}_all")
    feature_mat_arr=("features_rs,features_wm,features_lang,features_gamb,features_social,features_motor" \
            "features_cv,features_tbss_AD,features_schaefer_streamlen,features_lang" \
            "features_ct,features_tbss_OD,features_schaefer_streamcount_log,features_wm" \
            "features_cv,features_tbss_AD,features_schaefer_OD,features_lang" \
            "features_ct,features_cv,features_ca,features_rs,features_wm,features_lang, \
            features_gamb,features_social,features_motor,features_tbss_FA,features_tbss_MD, \
            features_tbss_AD,features_tbss_RD,features_tbss_OD,features_tbss_ICVF, \
            features_tbss_ISOVF,features_schaefer_streamlen,features_schaefer_FA, \
            features_schaefer_MD,features_schaefer_AD,features_schaefer_RD, \
            features_schaefer_OD,features_schaefer_ISOVF, \
            features_schaefer_ICVF,features_schaefer_streamcount_log")

## stats
elif [ $regression == "stats" ]; then
    base_cmd="$scriptdir/stats/CBIG_MMP_HCP_stats_perm_wrapper.sh "
    time_allowed="30:00:00"
    scheduler_tag="HCP_stats"
    mem="10GB"

    # additional parameters needed for stacking model
    stats_outdir=$projectdir/stats
    run_once=0
    # multikrr outstem
    multiKRR_outstem='multiKRR_k1_rs_k2_wm_k3_lang_k4_gamb_k5_social_k6_motor'
    # stacking outstem and features
    outstem_arr=("stacking_LRR_rs_wm_lang_gamb_social_motor" "stacking_LRR_all")
    feature_mat_arr=("features_rs,features_wm,features_lang,features_gamb,features_social,features_motor" \
            "features_ct,features_cv,features_ca,features_rs,features_wm,features_lang, \
            features_gamb,features_social,features_motor,features_tbss_FA,features_tbss_MD, \
            features_tbss_AD,features_tbss_RD,features_tbss_OD,features_tbss_ICVF, \
            features_tbss_ISOVF,features_schaefer_streamlen,features_schaefer_FA, \
            features_schaefer_MD,features_schaefer_AD,features_schaefer_RD, \
            features_schaefer_OD,features_schaefer_ISOVF, \
            features_schaefer_ICVF,features_schaefer_streamcount_log")
else
    echo "ERROR: type of regression not recognised, exiting..."
    exit
fi

# set up features to run
feature_array=("t1/features_ct" "t1/features_cv" "t1/features_ca" \
    "tbss/features_tbss_FA" "tbss/features_tbss_MD" "tbss/features_tbss_AD" \
    "tbss/features_tbss_RD" "tbss/features_tbss_OD" "tbss/features_tbss_ICVF" \
    "tbss/features_tbss_ISOVF" "tractography/features_schaefer_streamcount_log" \
    "tractography/features_schaefer_streamlen" "tractography/features_schaefer_FA" \
    "tractography/features_schaefer_MD" "tractography/features_schaefer_AD" \
    "tractography/features_schaefer_RD" "tractography/features_schaefer_OD" \
    "tractography/features_schaefer_ISOVF" "tractography/features_schaefer_ICVF" \
    "fmri/features_rs" "fmri/features_wm" "fmri/features_lang" "fmri/features_gamb" \
    "fmri/features_social" "fmri/features_motor") 

# common files and settings
outerFolds=10
innerFolds=10
subtxt=$behav_dir/MMP_HCP_753_subs.txt
# please modify path to score and restricted csv: NOT PROVIDED BY DEFAULT
scorecsv=$behav_dir/MMP_HCP_behavscores.csv
restrictedcsv=$behav_dir/RESTRICTED_jingweili_4_12_2017_1200subjects_fill_empty_zygosityGT_by_zygositySR.csv
predvar=$behav_dir/MMP_HCP_y_variables.txt
covtxt=$behav_dir/MMP_HCP_covariates.txt
ymat="MMP_HCP_753_y.mat"
covmat="MMP_HCP_753_covariates.mat"
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
        for split_idx in $(seq 1 1 60); do 
            cmd="$base_cmd $feature_path $outdir $outerFolds $innerFolds $subtxt $scorecsv \
            $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx "
            run_jobs
        done
    elif [ $regression == "LRR" ] || [ $regression == "elasticnet" ]; then
        feature_path=$feature_folder/$feature
        num_var=61
        for split_idx in $(seq 1 1 60); do 
            cmd="$base_cmd $feature_path $outdir $outerFolds $innerFolds $subtxt $scorecsv \
            $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx $num_var "
            run_jobs
        done
    elif [ $regression == "multiKRR" ]; then
        # run jobs to generate kernels for all splits first
        for split_idx in $(seq 1 1 60); do 
            cmd="$base_cmd $outstem $mkrr_feature_folder $outdir $outerFolds $innerFolds \
            $subtxt $scorecsv $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx 1 \
            $feature_mat_cells $kernel_groups "
            run_jobs
        done
        # remaining jobs only submitted if kernel for that split is already generated
        for split_idx in $(seq 1 1 60); do 
            for fold_idx in $(seq 2 1 10); do 
                # wait for kernels to be generated
                kernel_txt_file=$outdir/$outstem/seed_${split_idx}/results/step2_completion.txt
                while [ ! -f $kernel_txt_file ]; do
                    echo "Waiting for kernels to be generated..."
                    sleep 5m
                done
                cmd="$base_cmd $outstem $mkrr_feature_folder $outdir $outerFolds $innerFolds \
                    $subtxt $scorecsv $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx \
                    $fold_idx $feature_mat_cells $kernel_groups "
                run_jobs
            done
        done
        exit
    elif [ $regression == "stacking" ]; then
        num_var=61
        arraylength=${#outstem_arr[@]}
        for (( i=0; i<${arraylength}; i++ )); do
            outstem=${outstem_arr[$i]}
            feature_mat_cells=$( echo ${feature_mat_arr[$i]} | xargs )
            for split_idx in $(seq 1 1 60); do
                cmd="$base_cmd $outstem $krr_output_dir $outdir $outerFolds $innerFolds \
                $ymat $covmat $second_lvl $split_idx $num_var $feature_mat_cells "
                run_jobs
            done
        done
        exit
    elif [ $regression == "stats" ]; then
        for y_idx in $(seq 59 1 61); do
            # null models for single feature models
            featurebase=$(basename $feature)
            outstem="KRR_$featurebase"
            singlekrr_dir=$outdir/$outstem
            cmd="$base_cmd singleKRR $featurebase $stats_outdir $singlekrr_dir $subtxt $restrictedcsv $y_idx "
            run_jobs

            # null models for combined models
            if [[ $run_once == 0 ]]; then
                # null models for multiKRR models
                multikrr_dir=$outdir/$multiKRR_outstem
                cmd="$base_cmd multiKRR $multiKRR_outstem $stats_outdir $multikrr_dir $subtxt $restrictedcsv $y_idx "
                run_jobs

                # null models for stacking models
                arraylength=${#outstem_arr[@]}
                for (( i=0; i<${arraylength}; i++ )); do
                    stacking_outstem=${outstem_arr[$i]}
                    feature_mat_cells=$( echo ${feature_mat_arr[$i]} | xargs )
                    stacking_dir=$outdir/$stacking_outstem
                    cmd="$base_cmd stacking $stacking_outstem $stats_outdir $stacking_dir $subtxt \
                    $restrictedcsv $y_idx $outdir $feature_mat_cells "
                    run_jobs
                done
            fi
        done
        run_once=1
    fi
done
