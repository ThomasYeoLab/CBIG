#!/bin/bash
#
# This function runs the random forest regressions in ChenOoi2023.
# When calling this function, sample size should be specified as the first argument. 
# Sample size can be 400, 800, 2000, 3000 or 5260.
# The path to the folders containing the results should be specified as the second argument. 
#
# An example of how to call this function as follows:
#     ./CBIG_ICCW_submit_RF_train.sh 400 ~/storage/ICCW_replication
#
# Check the common variables section to see whether the data directories are correct.
# List of directories and files in common variables:
# 1. py_env: Name of python environment. Please ensure that the packages under `replication/config` are installed.
# 2. root_dir: The folder in which the random forest scripts are stored
# 3. out_dir: Output path for the results
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# sample size variable and results dir to be read in when running script
sample=$1
results_dir=$2

###################################################################################################
# set common variables
###################################################################################################
# directories (PLEASE MODIFY)
py_env="Chen2022_WR" 
root_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/RF
out_dir=$results_dir/RF

###################################################################################################
# run RF
###################################################################################################
for behav in $(seq 1 1 1); do # only run for 1 behavior, for full set replace with $(seq 1 1 39)
    for fold in $(seq 1 1 252); do
        # create file for output
        r_dir=$out_dir/$sample/behav_$behav/rng0/fold_$fold
        if [ ! -d $r_dir ]; then mkdir -p $r_dir; fi

        # do not submit job if .sav for model already exists
        if [ ! -f $r_dir/ABCD_RF_behav_$behav.sav ]; then
            log_f=$r_dir/RF_${sample}_${behav}_${fold}_train.log
            cmd="source activate $py_env; python ${root_dir}/CBIG_ICCW_RF_train.py \
                $fold $behav $sample $results_dir"
            ssh headnode "$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd '${cmd}' \
                    -walltime 8:00:00 -name 'ABCD_RF_train' \
                    -mem '30G' -joberr '$r_dir' -jobout '$log_f'" < /dev/null
        else
            echo "S:$sample, B:$behav, F:$fold --- Model exists!"
        fi
    done
done