#!/bin/sh
#####
# This script calls the matlab function to run LRR. User needs to provide the following variables.
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
# 11. y_idx: index of y variable to run regression for
# 
# EXAMPLE: 
#    CBIG_MMP_ABCD_LRR.sh $feature_path $outdir $sites $innerFolds $subtxt $subcsv \
#        $predvar $covtxt $ymat $covmat #y_idx
# EXAMPLE OF HOW TO CALL FUNCTION:
#    CBIG_MMP_ABCD_LRR.sh data_dir/features/rs.mat output_dir 3 10 data_dir/subs.txt \
#        data_dir/scores.csv data_dir/prediction_variables.txt \
#        data_dir/covariates.txt "output_y.mat" "output_cov.mat" 1
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

# set directories
script_dir=$(dirname "$(readlink -f "$0")")

# set outstem
feature_path=$1
featurebase=$(basename $feature_path)
outstem="LRR_$featurebase"

# set params
outdir=$2
sites=$3
innerFolds=$4
subtxt=$5
subcsv=$6
predvar=$7
covtxt=$8
ymat=$9
covmat=${10}
y_idx=${11}

# Create log file and save params
mkdir -p $outdir/$outstem/logs
LF="$outdir/$outstem/logs/${outstem}_${y_idx}.log"
if [ -f $LF ]; then rm $LF; fi
echo "outstem = $outstem" >> $LF
echo "outdir = $outdir" >> $LF
echo "sites = $sites" >> $LF
echo "innerFolds = $innerFolds" >> $LF
echo "subtxt = $subtxt" >> $LF
echo "subcsv = $subcsv" >> $LF
echo "predvar = $predvar" >> $LF
echo "covtxt = $covtxt" >> $LF
echo "ymat = $ymat" >> $LF
echo "covmat = $covmat" >> $LF
echo "y_idx = $y_idx" >> $LF

# Call matlab function
matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_ABCD_LRR( \
   $sites, $innerFolds, '$feature_path', '$featurebase', '$outdir', '$subtxt', '$subcsv', \
   '$predvar', '$covtxt', '$ymat', '$covmat', $y_idx); catch ME; display(ME.message); \
   end; exit; " >> $LF 2>&1
