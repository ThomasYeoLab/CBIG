#!/bin/sh
#####
# This function schedules the Haufe-inverted feature importance calculation in the CBIG scheduler. 
# EXAMPLE: 
#    CBIG_MMP_ABCD_interpretation_wrapper.sh $parentdir/replication
#
# First input is the project directory where results are being saved, for example:
# "/home/leon_ooi/storage/Multimodal_prediction_project/replication/ABCD"
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

### set up data directories
projectdir=$1
outdir=$projectdir/output
if [ -z $CBIG_ABCD_REPDATA_DIR ]; then
    echo "ABCD_REPDATA_DIR not found: using default location on CBIG HPC"
    export CBIG_ABCD_REPDATA_DIR=$CBIG_REPDATA_DIR/stable_projects/predict_phenotypes/Ooi2022_MMP/MMP_ABCD_data
fi
feature_folder=$CBIG_ABCD_REPDATA_DIR/data/features

logdir=$projectdir/scheduler_logs
if [ ! -d $logdir ]; then mkdir -p $logdir; fi
script_dir=$(dirname "$(readlink -f "$0")")

# set up features to run
feature_array=("t1/features_ct" "t1/features_cv" "t1/features_ca" \
    "tbss/features_tbss_FA" "tbss/features_tbss_MD" "tbss/features_tbss_AD" \
    "tbss/features_tbss_RD" "tbss/features_tbss_OD" "tbss/features_tbss_ICVF" \
    "tbss/features_tbss_ISOVF" "tractography/features_schaefer_streamcount_log" \
    "tractography/features_schaefer_streamlen" "tractography/features_schaefer_FA" \
    "tractography/features_schaefer_MD" "tractography/features_schaefer_AD" \
    "tractography/features_schaefer_RD" "tractography/features_schaefer_OD" \
    "tractography/features_schaefer_ISOVF" "tractography/features_schaefer_ICVF" \
    "fmri/features_rs" "fmri/features_mid" "fmri/features_sst" "fmri/features_nback")

### submit jobs to scheduler
for feature in ${feature_array[@]}; do
    feature_path=$feature_folder/$feature
    input_dir=$(dirname $feature_path ) 
    feature=$(basename $feature_path )
    cmd="$script_dir/CBIG_MMP_ABCD_Haufe.sh $input_dir $outdir $feature "
    ### define run job function: MODIFY THIS IF RUNNING ON YOUR OWN SCHEDULER
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime "24:00:00" -name "ABCD_Haufe" \
    -mem "12GB" -joberr "$logdir" -jobout "$logdir"
done