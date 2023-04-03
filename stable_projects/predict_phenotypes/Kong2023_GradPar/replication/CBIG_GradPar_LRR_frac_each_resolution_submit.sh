#!/bin/bash

# This script is used to submit jobs for each method and resolution for LRR fracridge.
# To save time, here we only submit jobs for Schaefer2018 in HCP dataset with 100 and 200 resolutions
# and for sICA in ABCD dataset with 50 and 100 resolutions. We will only run 3 splits for Schaefer2018.
# 
# Example:
#     sh CBIG_GradPar_LRR_frac_each_resolution_submit.sh ~/storage/Temporary/CBIG_GradPar_replication
#
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
out_dir=$1
if [ -d $out_dir ]; then
    rm -r $out_dir
fi
echo $out_dir
mkdir -p $out_dir/logs
out_dir=`realpath $1`
out_ABCD_dir=$out_dir/ABCD
out_HCP_dir=$out_dir/HCP

####################################################
# Submit jobs for each method and resolution for LRR
####################################################
# To save time, here we only submit jobs for Schaefer2018 in HCP dataset with 100 and 200 resolutions
# and for sICA in ABCD dataset with 50 and 100 resolutions. We will only run 3 splits for Schaefer2018
#
# To fully replicate Kong2023, please submit jobs for all methods and resolutions:
# Kong2021: 100,200,300,400,500,600,700,800,900,1000.
# Schaefer2018: 100,200,300,400,500,600,700,800,900,1000.
# PrincipalGrad: 1,5,10,20,40,60,80,100.
# LocalGrad: 1
# sICA: 50,100,200,300
# For HCP dataset, please run 100 splits.

# Schaefer2018 for HCP dataset
# Can set <curr_project_name> to Kong2021, Schaefer2018, PrincipalGrad, LocalGrad, sICA for different methods
curr_project_name='Schaefer2018'
# loop resolution 100,200:
# Can change the <res> range for each method
for (( res=100; res<=200; res+=100 )); do
    # Can change the <nsplit> range for each method
    for (( nsplit=1; nsplit<=3; nsplit+=1 )); do
        log_dir=$out_dir/logs/${curr_project_name}_${res}_${nsplit}.log

        cmd="cd ${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Kong2023_GradPar/replication;"
        cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
        cmd="$cmd \"clear;clc;close all;CBIG_GradPar_LRR_frac_each_resolution_HCP_wrapper \
        ${out_HCP_dir} ${curr_project_name} ${res} ${nsplit}; exit; \" "
        cmd="$cmd | tee -a $log_dir"

        ERR_FILE_PATH="${out_dir}/job_err_out/${curr_project_name}_${res}_${nsplit}/job.err"
        OUT_FILE_PATH="${out_dir}/job_err_out/${curr_project_name}_${res}_${nsplit}/job.out"
        mkdir -p ${out_dir}/job_err_out/${curr_project_name}_${res}_${nsplit}

        $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 3:00:00 -mem 16G \
        -name "GradPar_rep_${curr_project_name}_${res}_${nsplir}" -joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

        cmd=""
    done
done

# sICA for ABCD dataset
# Can set <curr_project_name> to Kong2021, Schaefer2018, PrincipalGrad, LocalGrad, sICA for different methods
curr_project_name='sICA'
# loop resolution 50, 100:
# Can change the <res> range for each method
for res in {50,100}; do
    log_dir=$out_dir/logs/${curr_project_name}_${res}.log

    cmd="cd ${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/Kong2023_GradPar/replication;"
    cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
    cmd="$cmd \"clear;clc;close all;CBIG_GradPar_LRR_frac_each_resolution_ABCD_wrapper \
    ${out_ABCD_dir} ${curr_project_name} ${res}; exit; \" "
    cmd="$cmd | tee -a $log_dir"

    ERR_FILE_PATH="${out_dir}/job_err_out/${curr_project_name}_${res}/job.err"
    OUT_FILE_PATH="${out_dir}/job_err_out/${curr_project_name}_${res}/job.out"
    mkdir -p ${out_dir}/job_err_out/${curr_project_name}_${res}

    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 3:00:00 -mem 16G \
    -name "GradPar_rep_${curr_project_name}_${res}" -joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

    cmd=""
done