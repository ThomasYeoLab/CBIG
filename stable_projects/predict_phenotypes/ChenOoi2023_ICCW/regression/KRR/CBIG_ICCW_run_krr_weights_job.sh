#!/bin/bash
# This function runs extracts regression weights for kernel regressions in ChenOoi2023.
# When calling this function, sample size should be specified as the first argument. 
# Sample size can be 400, 800, 2000, 3000 or 5260.
# The path to the folders containing the results should be specified as the second argument. 
#
# An example of how to call this function as follows:
#     ./CBIG_ICCW_run_krr_weights_job.sh 400 ~/storage/ICCW_replication
#
# Check the common variables section to see whether the data directories are correct.
# List of directories and files in common variables:
# 1. scripts_dir: The folder in which the KRR scripts are stored
# 2. krr_input: Features used for prediction (i.e. FC matrix for all subjects)
# 3. krr_result_dir: Directory in which KRR results were saved
# 4. outdir: output directory for results to be saved in
#
# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# sample size and results dir to be read in when running script
sample_size=$1
results_dir=$2

###################################################################################################
# set common variables
###################################################################################################
# directories (PLEASE MODIFY)
scripts_dir=${CBIG_CODE_DIR}/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/KRR
krr_input=${CBIG_REPDATA_DIR}/stable_projects/predict_phenotypes/ChenTam2022_TRBPC/FC/FC_subjects_rs_all_score_mf.mat
krr_result_dir=${results_dir}/KRR/${sample_size}
outdir=${results_dir}/KRR/${sample_size}/weights

###################################################################################################
# submit jobs to extract regression weights for each behavior
###################################################################################################
mkdir -p $outdir
for behav in `seq 1 1` # only run for 1 behavior, for full set replace with `seq 1 39`
do
    cmd="${scripts_dir}/CBIG_ICCW_KRR_weights.sh ${krr_input} ${krr_result_dir} ${behav} ${outdir}"
    errfile=${outdir}/job_err_out/weights_${behav}.err
    outfile=${outdir}/job_err_out/weights_${behav}.out
    ${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 2:00:0 -mem 32gb -joberr ${errfile} -jobout ${outfile}\
         -cmd "${cmd}" -name krr_weights
done
