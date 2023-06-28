#!/bin/bash

# Example:
#     sh CBIG_gwMRF_unit_test.sh ~/storage/Temporary/CBIG_gwMRF_unit_test
#
# Written by Yang Qing and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

user=`whoami`

###############################
# Set paths and create log file
###############################

# code_dir contains Schaefer2018 parcellation codes
code_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code

# data_dir contains data used by unit test
data_dir=$CBIG_TESTDATA_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/data

# ref_dir contains reference results
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


###################################
# Generate example results
###################################

ex_output_dir=${output_dir}/example_out
ex_input_fullpath=${ex_output_dir}/example_input_fullpaths.csv

# specify input parameters
start_idx='1'
end_idx='2'
num_left_cluster='50'
num_right_cluster='50'
smoothcost='5000'
num_iterations='2'
num_runs='1'
start_gamma='50000000'
exponential='15'

my_cmd="cd ${code_dir};"
my_cmd="${my_cmd} matlab -nosplash -nodisplay -nodesktop -r " 
my_cmd="${my_cmd} \"CBIG_gwMRF_build_data_and_perform_clustering ${ex_input_fullpath} ${ex_output_dir} \
 ${start_idx} ${end_idx} ${num_left_cluster} ${num_right_cluster} ${smoothcost} ${num_iterations} \
 ${num_runs} ${start_gamma} ${exponential}; exit;\" "

# submit job for generating example results
ERR_FILE_PATH="${ex_output_dir}/CBIG_gwMRF_example.err"
OUT_FILE_PATH="${ex_output_dir}/CBIG_gwMRF_example.out"

# this is to accomodate for new HPC
temp_script_file1="${ex_output_dir}/temp_script1.sh"
echo '#!/bin/bash ' >> ${temp_script_file1}
echo ${my_cmd} >> ${temp_script_file1}
chmod 755 ${temp_script_file1}

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script_file1}" -walltime 3:00:00 -mem 35G -name "gwMRF_ex" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

# will not proceed to next step until the current job is successfully done
ex_out_name='Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_5000_iterations_2_seed_1.mat'
ex_out_file=${ex_output_dir}/clustering/${ex_out_name}
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n gwMRF_ex -o ${ex_out_file}

EXIT_STATUS=$?
if [[ "$EXIT_STATUS" -eq "1" ]]
then
    echo "Failed to generate example parcellations."
    exit 1
else
    echo "Successfully generated example results. Will proceed with the rest of the unit test..."
fi

####################################
# Create input time and mult matrix
####################################

# note that the mult matrix is generated based on the time matrix
# therefore we squash two matlab commands together
# we only checks the final mult matrix result, no need to check time matrix (extra complexity)

TIME_DATA_OUT="${output_dir}/time_data"
MULT_DATA_OUT="${output_dir}/mult_mat"

mkdir -p ${TIME_DATA_OUT}
mkdir -p ${MULT_DATA_OUT}

# specifiy job inputs for build time mat
input_fullpaths="${output_dir}/unit_tests_input_fullpaths_short.csv"
start_idx='1'
end_idx='2'
fsaverage='fsaverage6'
lh_output_file='lh_time_matrix.mat'
rh_output_file='rh_time_matrix.mat'

# specify job inputs for build mult mat
lh_input_filename="${output_dir}/time_data/lh_time_matrix.mat"
rh_input_filename="${output_dir}/time_data/rh_time_matrix.mat"
lh_output_filename='lh_mult_matrix.mat'
rh_output_filename='rh_mult_matrix.mat'

# create one matlab command that sequencially calls the 2 functions
# that build time matrix and mult matrix separately
my_cmd2="cd ${code_dir}/lib ;"
my_cmd2="${my_cmd2} matlab -nosplash -nodisplay -nodesktop -r "
my_cmd2="${my_cmd2} \"CBIG_gwMRF_build_time_matrix ${input_fullpaths} ${TIME_DATA_OUT} ${start_idx} \
 ${end_idx} ${fsaverage} ${lh_output_file} ${rh_output_file}; "
my_cmd2="${my_cmd2} CBIG_gwMRF_build_prod_matrix ${lh_input_filename} ${rh_input_filename} \
${MULT_DATA_OUT} ${lh_output_filename} ${rh_output_filename}; exit; \" "

# submit job for generating the time matrix
ERR_FILE_PATH="${output_dir}/job_err_out/CBIG_gwMRF_mult.err"
OUT_FILE_PATH="${output_dir}/job_err_out/CBIG_gwMRF_mult.out"

# this is to accomodate for new HPC
temp_script_file2="${output_dir}/temp_script2.sh"
echo '#!/bin/bash' >> ${temp_script_file2}
echo ${my_cmd2} >> ${temp_script_file2}
chmod 755 ${temp_script_file2}

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script_file2}" -walltime 1:00:00 -mem 15G -name "gwMRF_mult_mat" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

# wait until the current job is done, by checking both the job name and the mult mat file
# we only check rh mult mat since it is generated after the lh mult mat
OUT_FILE_DIR="${output_dir}/mult_mat/rh_mult_matrix.mat"
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n 'gwMRF_mult_mat' -o ${OUT_FILE_DIR}

EXIT_STATUS=$?
if [[ "$EXIT_STATUS" -eq "1" ]]
then
    echo "[FAILED] Failed to generate multiplied matrices."
    exit 1
else
    echo "Successfully generated multiplied matrices. Will proceed with the rest of the unit test..."
fi

####################################################
# Generate intermediate RH parcellation for seed 2
####################################################

# in the function input below, note that we delibrately set 'iter_reduce_gamma' to 1
# since running the job to the end (iter_reduce_gamma = 300) will be very time consuming
# in the previous example, we have already run seed 1 to the end
# here, we only intend to test 1 more seed, thus 1 iteration is enough

vargin=" 'lh_avg_file' '${output_dir}/mult_mat/lh_mult_matrix.mat' \
  'rh_avg_file' '${output_dir}/mult_mat/rh_mult_matrix.mat' 'left_cluster' '50' 'right_cluster' \
  '50' 'iterations' '2' 'smoothcost' '1000' 'start_index' '2' 'runs' '1' 'output_folder' \
  '${output_dir}/clustering/' 'dim' '360' 'exponential' '15' 'start_gamma' '10000000' 'skip_left' '1' \
  'iter_reduce_gamma' '1' ";

PROJ_CODE_DIR="${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code"

my_cmd3="cd ${code_dir}; "
my_cmd3="${my_cmd3} matlab -nosplash -nodisplay -nodesktop -r "
my_cmd3="${my_cmd3} \" addpath(genpath('${PROJ_CODE_DIR}'));"
my_cmd3="${my_cmd3} CBIG_gwMRF_graph_cut_clustering ${vargin}; exit; \" ";

# submit job for generating the time matrix
ERR_FILE_PATH="${output_dir}/job_err_out/CBIG_gwMRF_seed2.err"
OUT_FILE_PATH="${output_dir}/job_err_out/CBIG_gwMRF_seed2.out"

# this is to accomodate for new HPC
temp_script_file3="${output_dir}/temp_script3.sh"
echo '#!/bin/bash' >> ${temp_script_file3}
echo ${my_cmd3} >> ${temp_script_file3}
chmod 755 ${temp_script_file3}
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script_file3}" -walltime 1:00:00 -mem 15G -name "gwMRF_seed2" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

##############################################
# Hold until first iteration results generated
##############################################
prefix=Graph_Cut_faster__grad_prior_gordon_cluster_20_datacost_1_smoothcost_1000_iterations_2
file_seed2=${output_dir}/clustering/inbetween_results/${prefix}_seed_2_rh_reduction_iteration_1.mat
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n 'gwMRF_seed2' -o ${file_seed2}

EXIT_STATUS=$?
if [[ "$EXIT_STATUS" -eq "1" ]]
then
    echo "[FAILED] Failed to generate results for seed 2."
    exit 1
else
    echo "Successfully generated results for seed 2. Will proceed with the rest of the unit test..."
fi

######################################################
# Compare results to check whether code pass unit test
######################################################

ERR_FILE_PATH=${output_dir}/job_err_out/CBIG_gwMRF_compare.err
OUT_FILE_PATH=${output_dir}/job_err_out/CBIG_gwMRF_compare.out

my_cmd4="cd `dirname $0`;"
my_cmd4="${my_cmd4} matlab -nosplash -nodisplay -nodesktop -r "
my_cmd4="${my_cmd4} \" clear;clc;close all;CBIG_gwMRF_check_unit_test_result '${output_dir}';exit; \" >> $LF "

# this is to accomodate for new HPC
temp_script_file4="${output_dir}/temp_script4.sh"

echo '#!/bin/bash' >> ${temp_script_file4}
echo ${my_cmd4} >> ${temp_script_file4}
chmod 755 ${temp_script_file4}
$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script_file4}" -walltime 1:00:00 -mem 30G -name "gwMRF_compare" \
-joberr ${ERR_FILE_PATH} -jobout ${OUT_FILE_PATH}

# hold until job is finished
sh $CBIG_CODE_DIR/utilities/scripts/CBIG_check_job_status -n 'gwMRF_compare'

# remove temporary scripts for pb_submit
rm ${temp_script_file1} 
rm ${temp_script_file2} 
rm ${temp_script_file3} 
rm ${temp_script_file4} 

# check if job contains no failure | is run to the end
fail=`grep FAILED $LF | wc -l`
done=`grep DONE $LF | wc -l`
if [ $fail -gt 0 ] || [ $done -eq 0 ]; then
    echo "[FAILED] Job is unsuccessful. Check log file: $LF."
    exit 1
fi

##################
# End of unit test
##################

echo -e "\nEnd CBIG_gwMRF_unit_test" >> $LF
date >> $LF
echo -e "====================================================================\n\n" >> $LF