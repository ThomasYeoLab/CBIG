#!/bin/bash

# This wrapper function performs KRR using a publicly-available air quality dataset.
# The aim of this example to let the users familiarize our KRR workflow, please
# refrain from making any conclusion regarding this air quality dataset per se.
# This script will firstly generate a setup file, then pass the setup file to KRR workflow.
# In this example, KRR is run using 5-fold cross validation with 3 different splits
#
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

outdir=$1

# prepare output directory
if [[ -d $outdir ]]
then
    rm -r $outdir
    mkdir $outdir
else
    mkdir $outdir
fi

for split in 1 2 3
do
    # specify path for KRR input
    input_path="${CBIG_CODE_DIR}/utilities/matlab/predictive_models/KernelRidgeRegression/example/input"
    
    # generate setup file for KRR
    cmd="matlab -nosplash -nodesktop -nodisplay -r \"CBIG_KRR_example_prepare_setup_file $input_path \
    $outdir $split; exit;\" "
    eval $cmd   
    setup_file_path="${outdir}/${split}/setup_${split}.mat"

    # pass the setup file to KRR workflow LITE
    cmd="matlab -nosplash -nodesktop -nodisplay -r \"CBIG_KRR_workflow_LITE $setup_file_path; exit;\" "
    eval $cmd 

done

# compare optimal accuracies
cmd="matlab -nosplash -nodesktop -nodisplay -r \"CBIG_KRR_example_check_result $outdir; exit;\" "
eval $cmd


exit 0
