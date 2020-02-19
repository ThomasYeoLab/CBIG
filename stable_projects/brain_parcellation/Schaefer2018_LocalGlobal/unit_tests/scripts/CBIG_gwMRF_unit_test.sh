#!/bin/bash

# Example:
#     sh CBIG_gwMRF_unit_test.sh ~/storage/Temporary/CBIG_gwMRF_unit_test
#
# Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

user=`whoami`


###############################
# Set paths and create log file
###############################

# code_dir containing Schaefer2018 parcellation codes
code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code

# data_dir containing data used by unit test
data_dir=$CBIG_TESTDATA_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/data

# ref_dir containing reference results
ref_dir=$CBIG_TESTDATA_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/results

# pass in output_dir
output_dir=`realpath ${1}`
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
# Create input time and mult matrix
####################################

$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ
#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'gwMRF_mult'
#PBS -l walltime=1:00:00
#PBS -l mem=15gb   
#PBS -e ${output_dir}/job_err_out/CBIG_gwMRF_mult.err
#PBS -o ${output_dir}/job_err_out/CBIG_gwMRF_mult.out

cd ${code_dir}/lib

mkdir -p ${output_dir}/time_data
mkdir -p ${output_dir}/mult_mat

matlab -nosplash -nodisplay -nodesktop -r "CBIG_gwMRF_build_time_matrix \
  '${data_dir}/unit_tests_input_fullpaths_short.csv' '${output_dir}/time_data/' '1' '2' \
  'fsaverage6' 'lh_time_matrix.mat' 'rh_time_matrix.mat'; \
  CBIG_gwMRF_build_prod_matrix '${output_dir}/time_data/lh_time_matrix.mat' \
  '${output_dir}/time_data/rh_time_matrix.mat' '${output_dir}/mult_mat/' \
  'lh_mult_matrix.mat' 'rh_mult_matrix.mat'; exit;" > /dev/null

EOJ

###########################################
# Hold until time and mult matrix generated
###########################################

flag=1
while [ "${flag}" == "1" ]
do 
	num_job=`qstat -u ${user} | grep gwMRF_mult | wc -l`
	if [ "${num_job}" == "1" ]; then
		sleep 3s
	else
		flag=0
	fi
done

job1_out_file=${output_dir}/mult_mat/rh_mult_matrix.mat
if [ ! -f $job1_out_file ]; then
	echo "[FAILED] Job is unsuccessful. Output file: $job1_out_file does not exist."
	exit 1
fi



####################################################
# Generate intermediate LH parcellation for seed 1
####################################################

$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ
#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'gwMRF_seed1'
#PBS -l walltime=1:00:00
#PBS -l mem=15gb   
#PBS -e ${output_dir}/job_err_out/CBIG_gwMRF_seed1.err
#PBS -o ${output_dir}/job_err_out/CBIG_gwMRF_seed1.out

cd ${code_dir}
matlab -nosplash -nodisplay -nodesktop -r "addpath(genpath(\
  '${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code')); \
  CBIG_gwMRF_graph_cut_clustering 'lh_avg_file' '${output_dir}/mult_mat/lh_mult_matrix.mat' \
  'rh_avg_file' '${output_dir}/mult_mat/rh_mult_matrix.mat' 'left_cluster' '50' 'right_cluster' \
  '50' 'iterations' '2' 'smoothcost' '1000' 'start_index' '1' 'runs' '1' 'output_folder' \
  '${output_dir}/clustering/' 'dim' '360' 'exponential' '15' 'start_gamma' '10000000' 'skip_right' '1';exit;"

EOJ

####################################################
# Generate intermediate RH parcellation for seed 2
####################################################

$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ
#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'gwMRF_seed2'
#PBS -l walltime=1:00:00
#PBS -l mem=15gb   
#PBS -e ${output_dir}/job_err_out/CBIG_gwMRF_seed2.err
#PBS -o ${output_dir}/job_err_out/CBIG_gwMRF_seed2.out

cd ${code_dir}
matlab -nosplash -nodisplay -nodesktop -r "addpath(genpath(\
  '${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code')); \
  CBIG_gwMRF_graph_cut_clustering 'lh_avg_file' '${output_dir}/mult_mat/lh_mult_matrix.mat' \
  'rh_avg_file' '${output_dir}/mult_mat/rh_mult_matrix.mat' 'left_cluster' '50' 'right_cluster' \
  '50' 'iterations' '2' 'smoothcost' '1000' 'start_index' '2' 'runs' '1' 'output_folder' \
  '${output_dir}/clustering/' 'dim' '360' 'exponential' '15' 'start_gamma' '10000000' 'skip_left' '1';exit;"

EOJ

##############################################
# Hold until first iteration results generated
##############################################
prefix=Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_2
file_seed1=${output_dir}/clustering/inbetween_results/${prefix}_seed_1_lh_reduction_iteration_1.mat
file_seed2=${output_dir}/clustering/inbetween_results/${prefix}_seed_2_rh_reduction_iteration_1.mat

flag=1
while [ "${flag}" == "1" ]
do 
	num_job=`qstat -u ${user} | grep gwMRF_seed1 | wc -l`
	if [ "${num_job}" == "1" ]; then
		
		if [ -f "${file_seed1}" ]; then
			job_id=`qstat -u ${user}| grep gwMRF_seed1 | awk '{print $1}' | sed 's/.circuv//g'`
			qdel ${job_id}
			flag=0
		else
			sleep 3s
		fi
		
	else
		flag=0
	fi
done

if [ ! -f $file_seed1 ]; then
	echo "[FAILED] Job is unsuccessful. Output file: $file_seed1 does not exist."
	exit 1
fi

flag=1
while [ "${flag}" == "1" ]
do 
	num_job=`qstat -u ${user} | grep gwMRF_seed2 | wc -l`
	if [ "${num_job}" == "1" ]; then
		
		if [ -f "${file_seed2}" ]; then
			job_id=`qstat -u ${user}| grep gwMRF_seed2 | awk '{print $1}' | sed 's/.circuv//g'`
			qdel ${job_id}
			flag=0
		else
			sleep 3s
		fi
		
	else
		flag=0
	fi
done

if [ ! -f $file_seed2 ]; then
	echo "[FAILED] Job is unsuccessful. Output file: $file_seed2 does not exist."
	exit 1
fi

######################################################
# Compare results to check whether code pass unit test
######################################################


$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ
#!/bin/sh
#PBS -S /bin/bash
#PBS -N 'gwMRF_compare'
#PBS -l walltime=1:00:00
#PBS -l mem=15gb   
#PBS -e ${output_dir}/job_err_out/CBIG_gwMRF_compare.err
#PBS -o ${output_dir}/job_err_out/CBIG_gwMRF_compare.out

cd `dirname $0` # cd to unit test scripts folder
matlab -nosplash -nodisplay -nodesktop -r "clear;clc;close all;CBIG_gwMRF_check_unit_test_result \
  '${output_dir}';exit" >> $LF
		
EOJ

flag=1
while [ "${flag}" == "1" ]
do 
	num_job=`qstat -u ${user} | grep gwMRF_compare | wc -l`
	if [ "${num_job}" == "1" ]; then
		sleep 3s
	else
		flag=0
	fi
done

fail=`grep FAILED $LF | wc -l`
success=`grep SUCCESS $LF | wc -l`
if [ $fail -gt 0 ] || [ $success -eq 0 ]; then
	echo "[FAILED] Job is unsuccessful. Check log file: $LF."
	exit 1
fi
		
##################
# End of unit test
##################

echo -e "\nEnd CBIG_gwMRF_unit_test" >> $LF
date >> $LF
echo -e "====================================================================\n\n" >> $LF








