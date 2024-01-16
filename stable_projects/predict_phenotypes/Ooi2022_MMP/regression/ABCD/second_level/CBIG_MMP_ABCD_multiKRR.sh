#!/bin/sh
#####
# This script calls the matlab function to run multiKRR. User needs to provide the following variables.
# 1. outstem: name of output file
# 2. mkrr_feature_folder: path to folder containing feature files
# 3. outdir: output directory
# 4. sites: number of sites used for the test fold
# 5. innerFold: number of inner folds
# 6. subtxt: list of subjects
# 7. subcsv: table of behaviour scores
# 8. predvar: txt file of names of behaviours to predict from subcsv
# 9. covtxt: txt file of names of covariates to regress from y variables
# 10. ymat: output name of behaviours to be predicted
# 11. covmat: output name of covariates to control for
# 12. fold_idx: index of fold to run KRR for
# 13. feature_mat_cells: string of comma separated file names for each feature
# 14. kernel_groups: string of comma separated groupings for multiKRR
# 
# EXAMPLE: 
#    CBIG_MMP_ABCD_multiKRR.sh $outstem $mkrr_feature_folder $outdir $sites $innerFolds \
#         $subtxt $subcsv $predvar $covtxt $ymat $covmat $fold_idx \
#         $feature_mat_cells $kernel_groups 
# EXAMPLE OF HOW TO CALL FUNCTION:
#    CBIG_MMP_ABCD_KRR.sh "multiKRR_rs_nback" data_dir 3 10 data_dir/subs.txt \
#        data_dir/scores.csv data_dir/prediction_variables.txt \
#        data_dir/covariates.txt "output_y.mat" "output_cov.mat" 1 \
#        "features_rs,features_nback" "1,2"
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

# set directories
script_dir=$(dirname "$(readlink -f "$0")")

# set outstem
outstem=$1

# set params
mkrr_feature_folder=$2
outdir=$3
sites=$4
innerFolds=$5
subtxt=$6
subcsv=$7
predvar=$8
covtxt=$9
ymat=${10}
covmat=${11}
fold_idx=${12}

# format feature array for matlab
# read features
feature_mat_cells="{"
IFS=',' read -ra delimited_list <<< ${13}
for i in "${delimited_list[@]}"; do
    feature_mat_cells=${feature_mat_cells}"'"$i"',"
done
feature_mat_cells=${feature_mat_cells::-1}"}"
# read kernel groups
kernel_groups="{"
IFS=',' read -ra delimited_list <<< ${14}
for i in "${delimited_list[@]}"; do
    kernel_groups=${kernel_groups}"["$i"],"
done
kernel_groups=${kernel_groups::-1}"}"


# Create log file and save params
mkdir -p $outdir/$outstem/logs
LF="$outdir/$outstem/logs/${outstem}_${fold_idx}.log"
if [ -f $LF ]; then rm $LF; fi
echo "outstem = $outstem" >> $LF
echo "mkrr_feature_folder = $mkrr_feature_folder" >> $LF
echo "outdir = $outdir" >> $LF
echo "sites = $sites" >> $LF
echo "innerFolds = $innerFolds" >> $LF
echo "subtxt = $subtxt" >> $LF
echo "subcsv = $subcsv" >> $LF
echo "predvar = $predvar" >> $LF
echo "covtxt = $covtxt" >> $LF
echo "ymat = $ymat" >> $LF
echo "covmat = $covmat" >> $LF
echo "fold_idx = $fold_idx" >> $LF
echo "feature_mat_cells = $feature_mat_cells" >> $LF
echo "kernel_groups = $kernel_groups" >> $LF

# Call matlab function
matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_ABCD_multiKRR($sites, \
    $innerFolds, '$mkrr_feature_folder', $feature_mat_cells, $kernel_groups, '$outstem', '$outdir', \
    '$subtxt', '$subcsv', '$predvar', '$covtxt', '$ymat', '$covmat', $fold_idx); catch ME; display(ME.message); \
    end; exit; " >> $LF 2>&1

