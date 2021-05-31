#!/bin/bash

# Example:
#     sh CBIG_ArealMSHBM_replication_wrapper.sh ~/storage/Temporary/CBIG_MSHBM_replication_wrapper
#
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
out_dir=`realpath $1`
if [ -d $out_dir ]; then
    rm -r $out_dir
fi
mkdir -p $out_dir

##################################################
# Generate functional connectivity profile list and 
# initialization parameters into output directory
##################################################

sh $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/replication\
/CBIG_ArealMSHBM_create_replication_input_data.sh $out_dir

#### Set up number of training subjects and number of iterations for the purpose of a simple replication
#### NOTE: To fully replicate Kong2022, set all_sub=40, num_iter="" 
all_sub=3
num_iter=2

#### Estimate group priors for HCP fs_LR_32k
code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step2_estimate_priors
project_dir_fslr=$out_dir/estimate_group_priors/HCP_fs_LR_32k
tmp_dir_fslr=$project_dir_fslr/tmp_results
mkdir -p $project_dir_fslr/logs

## Parent job
log_dir_fslr_parent=$project_dir_fslr/logs/HCP_fs_LR_32k_parent.log

cmd="cd ${code_dir};"
cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
cmd="$cmd \"clear;clc;close all;CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_parent \
${project_dir_fslr} fs_LR_32k $all_sub 4 400 100 ${tmp_dir_fslr} $num_iter; exit; \" "
cmd="$cmd | tee -a $log_dir_fslr_parent"

ERR_FILE_PATH="${project_dir_fslr}/job_err_out/parent/job.err"
OUT_FILE_PATH="${project_dir_fslr}/job_err_out/parent/job.out"
mkdir -p ${project_dir_fslr}/job_err_out/parent

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 72:00:00 -mem 72G -name "ArealMSHBM_rep_fslr_group_parent" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

cmd=""

## Children jobs
for (( sub=1; sub<=$all_sub; sub++ )); do
    log_dir_fslr_child=$project_dir_fslr/logs/HCP_fs_LR_32k_${sub}.log

    cmd="cd ${code_dir};"
    cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
    cmd="$cmd \"clear;clc;close all;CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_child \
    ${sub} fs_LR_32k ${tmp_dir_fslr}; exit; \" "
    cmd="$cmd | tee -a $log_dir_fslr_child"

    ERR_FILE_PATH="${project_dir_fslr}/job_err_out/child${sub}/job.err"
    OUT_FILE_PATH="${project_dir_fslr}/job_err_out/child${sub}/job.out"
    mkdir -p ${project_dir_fslr}/job_err_out/child${sub}

    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 72:00:00 -mem 8G \
    -name "ArealMSHBM_rep_fslr_group_child${sub}" -joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

    cmd=""
done

#### Estimate group priors for HCP fsaverage6
code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step2_estimate_priors
project_dir_fs6=$out_dir/estimate_group_priors/HCP_fsaverage6
tmp_dir_fs6=$project_dir_fs6/tmp_results
mkdir -p $project_dir_fs6/logs

## Parent job
log_dir_fs6_parent=$project_dir_fs6/logs/HCP_fsaverage6_parent.log

cmd="cd ${code_dir};"
cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
cmd="$cmd \"clear;clc;close all;CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_parent \
${project_dir_fs6} fsaverage6 ${all_sub} 4 400 5 ${tmp_dir_fs6} $num_iter; exit; \" "
cmd="$cmd | tee -a $log_dir_fs6_parent"

ERR_FILE_PATH="${project_dir_fs6}/job_err_out/parent/job.err"
OUT_FILE_PATH="${project_dir_fs6}/job_err_out/parent/job.out"
mkdir -p ${project_dir_fs6}/job_err_out/parent

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 72:00:00 -mem 72G -name "ArealMSHBM_rep_fs6_group_parent" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

cmd=""

## Children jobs
for (( sub=1; sub<=$all_sub; sub++ )); do
    log_dir_fs6_child=$project_dir_fs6/logs/HCP_fsaverage6_${sub}.log

    cmd="cd ${code_dir};"
    cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
    cmd="$cmd \"clear;clc;close all;CBIG_ArealMSHBM_gMSHBM_estimate_group_priors_child \
    ${sub} fsaverage6 ${tmp_dir_fs6}; exit; \" "
    cmd="$cmd | tee -a $log_dir_fs6_child"

    ERR_FILE_PATH="${project_dir_fs6}/job_err_out/child${sub}/job.err"
    OUT_FILE_PATH="${project_dir_fs6}/job_err_out/child${sub}/job.out"
    mkdir -p ${project_dir_fs6}/job_err_out/child${sub}

    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 72:00:00 -mem 8G \
    -name "ArealMSHBM_rep_fs6_group_child${sub}" -joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

    cmd=""
done


#### Generate individual parcellations for HCP fs_LR_32k
code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step3_generate_ind_parcellations
project_dir_fslr=$out_dir/generate_individual_parcellations/HCP_fs_LR_32k
log_dir_fslr=$project_dir_fslr/HCP_fs_LR_32k.log


cmd="cd ${code_dir};"
cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
cmd="$cmd \"clear;clc;close all;CBIG_ArealMSHBM_gMSHBM_generate_individual_parcellation \
${project_dir_fslr} fs_LR_32k 3 400 1 30 30 50 test_set; exit; \" "
cmd="$cmd | tee -a $log_dir_fslr"

ERR_FILE_PATH="${project_dir_fslr}/job_err_out/job.err"
OUT_FILE_PATH="${project_dir_fslr}/job_err_out/job.out"
mkdir -p ${project_dir_fslr}/job_err_out

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 01:00:00 -mem 6G -name "ArealMSHBM_rep_fslr_indpar" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

cmd=""

#### Generate individual parcellations for ABCD fsaverage6
code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2022_ArealMSHBM/step3_generate_ind_parcellations
project_dir_fs6=$out_dir/generate_individual_parcellations/ABCD_fsaverage6
log_dir_fs6=$project_dir_fs6/ABCD_fsaverage6.log

cmd="cd ${code_dir};"
cmd="$cmd matlab -nosplash -nodisplay -nodesktop -r "
cmd="$cmd \"clear;clc;close all;CBIG_ArealMSHBM_gMSHBM_generate_individual_parcellation \
${project_dir_fs6} fsaverage6 2 400 1 50 50 5 test_set; exit; \" "
cmd="$cmd | tee -a $log_dir_fs6"

ERR_FILE_PATH="${project_dir_fs6}/job_err_out/job.err"
OUT_FILE_PATH="${project_dir_fs6}/job_err_out/job.out"
mkdir -p ${project_dir_fs6}/job_err_out

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${cmd}" -walltime 01:00:00 -mem 6G -name "ArealMSHBM_rep_fs6_indpar" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}


#########################################################################
# Check if the parent and child jobs for estimating group priors crashed
#########################################################################
#
flag_fslr=1
flag_fs6=1
user_name=`(whoami)`
while [ "$flag_fslr" == "1" ] || [ "$flag_fs6" == "1" ]; do
    log_done_fslr=`(grep -l "Finished" $out_dir/estimate_group_priors/HCP_fs_LR_32k/logs/*)`
    job_err_fslr=`(find $out_dir/estimate_group_priors/HCP_fs_LR_32k/job_err_out -type f -name "job.err")`

    log_done_fs6=`(grep -l "Finished" $out_dir/estimate_group_priors/HCP_fsaverage6/logs/*)`
    job_err_fs6=`(find $out_dir/estimate_group_priors/HCP_fsaverage6/job_err_out -type f -name "job.err")`

    # if there is a job.err file in the job_err_out folder, at least one job is crashed/finished
    if [ -n "$job_err_fslr" ]; then
        # if at least one job is terminated, but the log file does not contain "Finished"
        # this suggests the job is crashed with error
        if [ -z "$log_done_fslr" ]; then
            echo "Please check your estimate_group_priors/HCP_fs_LR_32k/job_err_out, some jobs crashed!"
            echo "Kill the current running jobs ..."
            qdel $(qselect -u $user_name)
            exit 1
        else
            echo "Job finished without error!"
            flag_fslr=0
        fi
    fi

    # if there is a job.err file in the job_err_out folder, at least one job is crashed/finished
    if [ -n "$job_err_fs6" ]; then
        # if at least one job is terminated, but the log file does not contain "Finished"
        # this suggests the job is crashed with error
        if [ -z "$log_done_fs6" ]; then
            echo "Please check your estimate_group_priors/HCP_fsaverage6/job_err_out, some jobs crashed!"
            echo "Kill the current running jobs ..."
            qdel $(qselect -u $user_name)
            exit 1
        else
            echo "Job finished without error!"
            flag_fs6=1
        fi
    fi
done