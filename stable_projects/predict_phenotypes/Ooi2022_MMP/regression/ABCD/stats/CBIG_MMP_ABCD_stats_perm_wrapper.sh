#!/bin/sh
#####
# This script calls the matlab function to run generate the null models. User needs to provide the following variables.
# 1. feature_path: path to feature mat file
# 2. outdir: output directory
# 3. sites: number of sites used for the test fold
# 4. innerFold: number of inner folds
# 5. subtxt: list of subjects
# 6. subcsv: table of behaviour scores
# 7. predvar: txt file of names of behaviours to predict from subcsv
# 8. covtxt: txt file of names of covariates to regress from y variables
# 9. ymat: output name of behaviours to be predicted
# 10. covmat: output name of covariates to control for
# 
# EXAMPLE: 
#    CBIG_MMP_ABCD_KRR.sh $feature_path $outdir $sites $innerFolds \
#        $subtxt $subcsv $predvar $covtxt $ymat $covmat
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

# set directories
script_dir=$(dirname "$(readlink -f "$0")")

# set outstem
model=$1
# append outstem
if [ $model == "singleKRR" ]; then
    outstem="KRR_$2"
elif [ $model == "multiKRR" ]; then
    outstem=$2
elif [ $model == "stacking" ]; then
    outstem=$2
    singlekrr_outdir=$7
    feature_mat_cells="{"
    IFS=',' read -ra delimited_list <<< $8
    for i in "${delimited_list[@]}"; do
        feature_mat_cells=${feature_mat_cells}"'"$i"',"
    done
    feature_mat_cells=${feature_mat_cells::-1}"}"
fi

# set params
stats_outdir=$3
outdir=$stats_outdir/$outstem
model_dir=$4
site_list=$5
score_ind=$6
perm_seed_start=1
N_perm=10000

# Create log file and save params
mkdir -p $outdir/logs
LF="$outdir/logs/${outstem}_${score_ind}.log"
if [ -f $LF ]; then rm $LF; fi
echo "model = $model" >> $LF
echo "outstem = $outstem" >> $LF
echo "stats_outdir = $outdir" >> $LF
echo "model_dir = $model_dir" >> $LF
echo "site_list = $site_list" >> $LF
echo "score_ind = $score_ind" >> $LF
echo "perm_seed_start = $perm_seed_start" >> $LF
echo "N_perm = $N_perm" >> $LF

# Call matlab function
if [ $model == "singleKRR" ]; then
    matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_ABCD_compute_singleKRR_perm_stats(\
        '$model_dir', '$outstem', $score_ind, $perm_seed_start, $N_perm, '$site_list', \
        '$outdir'); catch ME; display(ME.message); end; exit; " >> $LF 2>&1

elif [ $model == "multiKRR" ]; then
    matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_ABCD_compute_multiKRR_perm_stats(\
        '$model_dir', '$outstem', $score_ind, $perm_seed_start, $N_perm, '$site_list', \
        '$outdir'); catch ME; display(ME.message); end; exit; " >> $LF 2>&1

elif [ $model == "stacking" ]; then
    matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_ABCD_compute_stacking_perm_stats(\
        '$singlekrr_outdir', '$model_dir', '$outstem', $feature_mat_cells, $score_ind, $perm_seed_start, $N_perm, \
        '$site_list', '$outdir'); catch ME; display(ME.message); end; exit; " >> $LF 2>&1

fi