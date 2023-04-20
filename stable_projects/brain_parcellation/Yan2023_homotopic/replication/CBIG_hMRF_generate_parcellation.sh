#!/bin/bash
# This function submits the job for generation of the final parcellation.
#
# Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################
# Set up output folder
###########################
out_dir=`realpath $1`
log_dir="${out_dir}/log_files_parc_generation"
mkdir -p ${log_dir}

# specify log/std out directories.
MANUAL_LOG_FILE="${log_dir}/manual_log.txt"
STD_OUT_FILE="${log_dir}/stdout.txt"
STD_ERR_FILE="${log_dir}/stderr.txt"

########################################
# Resolution specific input parameters #
########################################
start_seed=$2
num_cluster=$3
c=$4
d=$5
w_xyz=$6
decrease_d=$7
decrease_c=$8
premature_stopping=$9
lhrh_avg_file=${10}
job_name=${11}

#######################################
# Fixed Parameters across resolutions #
#######################################
k=15
seeds_per_job=1
num_iterations=100
mesh_type='fsaverage6'

#########################
# Submit job to cluster #
#########################
CODE_DIR="${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yan2023_homotopic/code/step2_generate_parcellation"
my_cmd="cd ${CODE_DIR};"
my_cmd="${my_cmd} matlab -nosplash -nodisplay -nodesktop -r \"clear;clc;close all;" 
my_cmd="${my_cmd} CBIG_hMRF_wrapper_generate_homotopic_parcellation ${lhrh_avg_file} ${out_dir} \
'start_seed' ${start_seed} 'num_rand_inits' ${seeds_per_job} 'num_cluster' ${num_cluster} 'num_iterations'\
 ${num_iterations} 'initial_c' ${c} 'initial_d' ${d} 'k' ${k} 'w_xyz' ${w_xyz} 'mesh_type' ${mesh_type} \
 'decrease_c' ${decrease_c} 'decrease_d' ${decrease_d} 'premature_stopping' ${premature_stopping}\
  ; exit;\" |tee -a ${MANUAL_LOG_FILE}"

# circumvent the bug ssh from compiler to head node
temp_script="${out_dir}/temp/temp_script_generate_parcellation.sh"
if [ -f "${temp_script}" ]; then
    rm ${temp_script}
fi
echo '#!/bin/bash' >> ${temp_script}
echo ${my_cmd} >> ${temp_script}
chmod 755 ${temp_script}

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script}" -walltime "30:00:00" -mem '75G' -name "${job_name}" \
-joberr ${STD_ERR_FILE} -jobout ${STD_OUT_FILE} 
# Do not remove double quotes. Else the interpreter cannot recognize my_cmd as an entire string