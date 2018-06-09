#!/bin/bash

# Example:
#     sh CBIG_gwMRF_unit_test.sh ~/storage/Temporary/CBIG_gwMRF_unit_test
#
# Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


###############################
# Set paths and create log file
###############################

# code_dir containing Schaefer2018 parcellation codes
code_dir=~/storage/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code

# data_dir containing data used by unit test
data_dir=/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/data

# ref_dir containing reference results
ref_dir=/mnt/eql/yeo1/CBIG_private_unit_tests_data/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/results

# pass in output_dir
output_dir=${1}
mkdir -p ${output_dir}/logs
mkdir -p ${output_dir}/job_err_out

# create log file in output_dir
LF="${output_dir}/logs/CBIG_gwMRF_unit_test.log"
touch $LF
echo -e "Start CBIG_gwMRF_unit_test" >> $LF
date >> $LF
echo -e "\nOutput folder = ${output_dir}" >> $LF
echo -e "Reference folder = ${ref_dir} \n\n" >> $LF


####################################
# Submit unit test job to circ-spool
####################################

/apps/sysapps/TORQUE/bin/qsub -V -q circ-spool << EOJ

#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'CBIG_gwMRF_unit_test'
#PBS -l walltime=15:00:00
#PBS -l mem=15gb   
#PBS -e ${output_dir}/job_err_out/CBIG_gwMRF_unit_test.err
#PBS -o ${output_dir}/job_err_out/CBIG_gwMRF_unit_test.out

####################################################
# Run unit test example and get parcellation results
####################################################

cd ${code_dir}
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_gwMRF_build_data_and_perform_clustering '${data_dir}/unit_tests_input_fullpaths_modified.csv' '${output_dir}' '1' '10' '50' '50' '1000' '7' '2' '10000000' '15';exit;"

######################################################
# Compare results to check whether code pass unit test
######################################################

cd `dirname $0` # cd to unit test scripts folder
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_gwMRF_check_unit_test_result '${output_dir}';exit" >> $LF

##################
# End of unit test
##################

echo -e "\nEnd CBIG_gwMRF_unit_test" >> $LF
date >> $LF
echo -e "====================================================================\n\n" >> $LF

EOJ






