#!/bin/sh
#####
# This script calls the matlab function to run multiKRR. User needs to provide the following variables.
# 1. outstem: name of output file
# 2. outdir: output directory
# 3. krr_output_dir: path to folder containing first level regression results
# 4. outerFolds: number of outer folds
# 5. innerFold: number of inner folds
# 6. ymat: output name of behaviours to be predicted
# 7. covmat: output name of covariates to control for
# 8. second_lvl: type of regression for stacking
# 9. split_idx: index of split
# 10. num_var: number of behaviour variables
# 11. feature_mat_cells: string of comma separated file names for each feature
# 
# EXAMPLE: 
#    CBIG_MMP_HCP_stacking.sh $outstem $krr_output_dir $outdir $outerFolds $innerFolds \
#        $ymat $covmat $second_lvl $split_idx $num_var $feature_mat_cells
# EXAMPLE OF HOW TO CALL FUNCTION:
#    CBIG_MMP_HCP_stacking.sh "stacking_LRR" output_dir output_dir 10 10 \
#        "output_y.mat" "output_cov.mat" LRR 1 61 "features_rs,features_wm"
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# set directories
script_dir=$(dirname "$(readlink -f "$0")")

# set outstem
outstem=$1

# set params
krr_output_dir=$2
outdir=$3
outerFolds=$4
innerFolds=$5
ymat=$6
covmat=$7
second_lvl=$8
split_idx=$9
num_var=${10}

# format feature array for matlab
# read features
# format for matlab
feature_mat_cells="{"
IFS=',' read -ra delimited_list <<< ${11}
for i in "${delimited_list[@]}"; do
    feature_mat_cells=${feature_mat_cells}"'"$i"',"
done
feature_mat_cells=${feature_mat_cells::-1}"}"

# loop over y variables
for y_idx in $(seq 1 1 $num_var); do
    # Create log file and save params
    mkdir -p $outdir/$outstem/logs
    LF="$outdir/$outstem/logs/${outstem}_${split_idx}_${y_idx}.log"
    if [ -f $LF ]; then rm $LF; fi
    echo "outstem = $outstem" >> $LF
    echo "krr_output_dir = $krr_output_dir" >> $LF
    echo "outdir = $outdir" >> $LF
    echo "outerFolds = $outerFolds" >> $LF
    echo "innerFolds = $innerFolds" >> $LF
    echo "ymat = $ymat" >> $LF
    echo "covmat = $covmat" >> $LF
    echo "second_lvl = $second_lvl" >> $LF
    echo "split_idx = $split_idx" >> $LF
    echo "y_idx = $y_idx" >> $LF
    echo "feature_mat_cells = $feature_mat_cells" >> $LF

    # Call matlab function
    matlab -nodesktop -nosplash -nodisplay -r " try addpath('$script_dir'); CBIG_MMP_HCP_stacking( \
        $outerFolds, $innerFolds, '$krr_output_dir', $feature_mat_cells, '$outstem', '$outdir', \
        $split_idx, '$ymat', '$covmat', $y_idx, '$second_lvl'); catch ME; display(ME.message); \
        end; exit; " >> $LF 2>&1
done
