#!/bin/bash
# This function submit the job for generation of premultiplied matrix.
#
# Input arguments:
# 1) out_dir: output directory to which the output matrices would be saved.
# 2) job_name: name of the job to be submitted to the cluster.
# 3) lh_fmri_fullpath_txt: left hemisphere subject fullpath input.
# 4) rh_fmri_fullpath_txt: right hemisphere subject fullpath input.
#
# Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

########################
# Set up output folder #
########################
out_dir=`realpath $1`
job_name=$2
lh_fmri_fullpath_txt=$3
rh_fmri_fullpath_txt=$4

log_dir="${out_dir}/log_file_PMM"
mkdir -p ${log_dir}

# specify log/std out directories.
MANUAL_LOG_FILE="${log_dir}/manual_log.txt"
STD_OUT_FILE="${log_dir}/stdout.txt"
STD_ERR_FILE="${log_dir}/stderr.txt"

#######################
# Specify input paths #
#######################
mesh_type='fsaverage6'

CODE_DIR="${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yan2023_homotopic/code/step1_generate_fmri_input"
my_cmd="cd ${CODE_DIR};"
my_cmd="${my_cmd} matlab -nosplash -nodisplay -nodesktop -r \"clear;clc;close all;" 
my_cmd="${my_cmd} CBIG_hMRF_generate_premultiplied_matrix ${out_dir} ${lh_fmri_fullpath_txt}\
 ${rh_fmri_fullpath_txt} ${mesh_type}; exit;\" |tee -a ${MANUAL_LOG_FILE}"

# circumvent the bug ssh from compiler to head node
temp_script="${out_dir}/temp/temp_script_generate_PMM.sh"
if [ -f "${temp_script}" ]; then
    rm ${temp_script}
fi
echo '#!/bin/bash' >> ${temp_script}
echo ${my_cmd} >> ${temp_script}
chmod 755 ${temp_script}

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script}" -walltime "10:00:00" -mem 200G -name ${job_name} \
-joberr ${STD_ERR_FILE} -jobout ${STD_OUT_FILE}
