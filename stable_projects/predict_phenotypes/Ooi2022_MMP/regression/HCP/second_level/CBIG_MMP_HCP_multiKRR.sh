#!/bin/sh
#####
# This script calls the matlab function to run multiKRR. User needs to provide the following variables.
# 1. outstem: name of output file
# 2. mkrr_feature_folder: path to folder containing feature files
# 3. outdir: output directory
# 4. outerFolds: number of outer folds
# 5. innerFold: number of inner folds
# 6. subtxt: list of subjects
# 7. scorecsv: table of behaviour scores and gender
# 8. restrictedcsv: table of family ID and age
# 9. predvar: txt file of names of behaviours to predict from subcsv
# 10. covtxt: txt file of names of covariates to regress from y variables
# 11. ymat: output name of behaviours to be predicted
# 12. covmat: output name of covariates to control for
# 13. split_idx: current split
# 14. fold_idx: current fold
# 15. feature_mat_cells: string of comma separated file names for each feature
# 16. kernel_groups: string of comma separated groupings for multiKRR
# 
# EXAMPLE: 
#    CBIG_MMP_HCP_multiKRR.sh $feature_path $outdir $outerFolds $innerFolds $subtxt $scorecsv \
#        $restrictedcsv $predvar $covtxt $ymat $covmat $split_idx
# EXAMPLE OF HOW TO CALL FUNCTION:
#    CBIG_MMP_HCP_KRR.sh "multiKRR_rs_wm" data_dir 10 10 data_dir/subs.txt \
#        data_dir/scores.csv data_dir/HCP_restricted.csv data_dir/prediction_variables.txt \
#        data_dir/covariates.txt "output_y.mat" "output_cov.mat" 1 1 \
#        "features_rs,features_wm" "1,2"
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# set directories
script_dir=$(dirname "$(readlink -f "$0")")

# set outstem
outstem=$1

# set params
mkrr_feature_folder=$2
outdir=$3
outerFolds=$4
innerFolds=$5
subtxt=$6
scorecsv=$7
restrictedcsv=$8
predvar=$9
covtxt=${10}
ymat=${11}
covmat=${12}
split_idx=${13}
fold_idx=${14}

# format feature array for matlab
# read features
feature_mat_cells="{"
IFS=',' read -ra delimited_list <<< ${15}
for i in "${delimited_list[@]}"; do
    feature_mat_cells=${feature_mat_cells}"'"$i"',"
done
feature_mat_cells=${feature_mat_cells::-1}"}"
# read kernel groups
kernel_groups="{"
IFS=',' read -ra delimited_list <<< ${16}
for i in "${delimited_list[@]}"; do
    kernel_groups=${kernel_groups}"["$i"],"
done
kernel_groups=${kernel_groups::-1}"}"

# Create log file and save params
mkdir -p $outdir/$outstem/logs
LF="$outdir/$outstem/logs/${outstem}_${split_idx}_fold${fold_idx}.log"
if [ -f $LF ]; then rm $LF; fi
echo "outstem = $outstem" >> $LF
echo "mkrr_feature_folder = $mkrr_feature_folder" >> $LF
echo "outdir = $outdir" >> $LF
echo "outerFolds = $outerFolds" >> $LF
echo "innerFolds = $innerFolds" >> $LF
echo "subtxt = $subtxt" >> $LF
echo "scorecsv = $scorecsv" >> $LF
echo "restrictedcsv = $restrictedcsv" >> $LF
echo "predvar = $predvar" >> $LF
echo "covtxt = $covtxt" >> $LF
echo "ymat = $ymat" >> $LF
echo "covmat = $covmat" >> $LF
echo "split_idx = $split_idx" >> $LF
echo "fold_idx = $fold_idx" >> $LF
echo "feature_mat_cells = $feature_mat_cells" >> $LF
echo "kernel_groups = $kernel_groups" >> $LF

# Call matlab function
matlab -nodesktop -nosplash -nodisplay -r "try addpath('$script_dir'); CBIG_MMP_HCP_multiKRR( \
    $outerFolds, $innerFolds, '$mkrr_feature_folder', $feature_mat_cells, $kernel_groups, \
    '$outstem', $split_idx, '$outdir', '$subtxt', '$scorecsv', '$restrictedcsv', '$predvar', \
    '$covtxt', '$ymat', '$covmat', $fold_idx); catch ME; display(ME.message); end; exit; " >> $LF 2>&1

