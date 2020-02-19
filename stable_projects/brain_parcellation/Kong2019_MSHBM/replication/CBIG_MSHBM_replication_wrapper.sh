#!/bin/bash

# Example:
#     sh CBIG_MSHBM_unit_test.sh ~/storage/Temporary/CBIG_MSHBM_replication_wrapper
#
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
out_dir=`realpath $1`
mkdir -p $out_dir

##################################################
# Generate functional connectivity profile list and 
# initialization parameters into output directory
##################################################

sh $CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/replication\
/CBIG_MSHBM_create_replication_input_data.sh $out_dir

code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Kong2019_MSHBM/step2_estimate_priors


project_dir_fs5=$out_dir/estimate_group_priors/GSP_HNU
project_dir_fslr=$out_dir/estimate_group_priors/HCP


log_dir_fs5=$project_dir_fs5/GSP.log
log_dir_fslr=$project_dir_fslr/HCP.log


###########################################
# Submit replication job to circ-spool
###########################################

# GSP
$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ

#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'CBIG_MSHBM_estimate_group_priors'
#PBS -l walltime=72:00:00
#PBS -l mem=120gb
#PBS -l nodes=4:ppn=4   
#PBS -e ${out_dir}/estimate_group_priors_GSP.err
#PBS -o ${out_dir}/estimate_group_priors_GSP.out


####################################
# unit test 1: estimate group priors
####################################


cd ${code_dir}
matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all;CBIG_MSHBM_estimate_group_priors '${project_dir_fs5}' 'fsaverage5' '37' '2' '17';exit;" \
>> $log_dir_fs5


echo -e "\nEnd CBIG_MSHBM_estimate_group_priors" >> $log_dir_fs5
date >> $log_dir_fs5
echo -e "====================================================================\n\n" >> $log_dir_fs5

EOJ

# HCP
$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ

#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'CBIG_MSHBM_estimate_group_priors'
#PBS -l walltime=72:00:00
#PBS -l mem=120gb
#PBS -l nodes=4:ppn=4   
#PBS -e ${out_dir}/estimate_group_priors_HCP.err
#PBS -o ${out_dir}/estimate_group_priors_HCP.out


####################################
# unit test 1: estimate group priors
####################################


cd ${code_dir}
matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all;CBIG_MSHBM_estimate_group_priors '${project_dir_fslr}' 'fs_LR_32k' '40' '4' '17';exit;" \
>> $log_dir_fslr


echo -e "\nEnd CBIG_MSHBM_estimate_group_priors" >> $log_dir_fslr
date >> $log_dir_fslr
echo -e "====================================================================\n\n" >> $log_dir_fslr

EOJ





