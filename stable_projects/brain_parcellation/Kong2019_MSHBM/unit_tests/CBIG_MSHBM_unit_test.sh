#!/bin/bash

# Example:
#     sh CBIG_MSHBM_unit_test.sh ~/storage/Temporary/CBIG_MSHBM_unit_test
#
# Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Specify output directory
###########################
out_dir=$1
mkdir -p $out_dir

##################################################
# Copy functional connectivity profile list and 
# initialization parameters into output directory
##################################################

# copy files for group priors estimation
cp -r /mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Kong2019_MSHBM/estimate_group_priors $out_dir

# copy files for individual parcellations generation
cp -r /mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Kong2019_MSHBM/generate_individual_parcellations $out_dir

code_dir1=/data/users/rkong/storage/ruby/data/Kong2019_MSHBM/step2_estimate_priors
code_dir2=/data/users/rkong/storage/ruby/data/Kong2019_MSHBM/step3_generate_ind_parcellations

project_dir1_fs5=$out_dir/estimate_group_priors/GSP_HNU
project_dir1_fslr=$out_dir/estimate_group_priors/HCP
project_dir2_fs5=$out_dir/generate_individual_parcellations/GSP_HNU
project_dir2_fslr=$out_dir/generate_individual_parcellations/HCP

log_dir1_fs5=$project_dir1_fs5/GSP.log
log_dir1_fslr=$project_dir1_fslr/HCP.log


###########################################
# Submit unit test job to circ-spool
###########################################

# GSP
/apps/sysapps/TORQUE/bin/qsub -V -q circ-spool << EOJ

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


cd ${code_dir1}
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_MSHBM_estimate_group_priors '${project_dir1_fs5}' 'fsaverage5' '37' '2' '17';exit;" >> $log_dir1_fs5


echo -e "\nEnd CBIG_MSHBM_estimate_group_priors" >> $log_dir1_fs5
date >> $log_dir1_fs5
echo -e "====================================================================\n\n" >> $log_dir1_fs5

EOJ

# HCP
/apps/sysapps/TORQUE/bin/qsub -V -q circ-spool << EOJ

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


cd ${code_dir1}
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_MSHBM_estimate_group_priors '${project_dir1_fslr}' 'fs_LR_32k' '40' '4' '17';exit;" >> $log_dir1_fslr


echo -e "\nEnd CBIG_MSHBM_estimate_group_priors" >> $log_dir1_fslr
date >> $log_dir1_fslr
echo -e "====================================================================\n\n" >> $log_dir1_fslr

EOJ





