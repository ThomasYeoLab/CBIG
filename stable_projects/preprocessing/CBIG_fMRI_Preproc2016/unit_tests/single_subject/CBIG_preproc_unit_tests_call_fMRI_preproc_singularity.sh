#!/bin/bash
# This script submit a job to HPC for CBIG_fMRI_Preproc2016 single subject unit test within a Singularity container.
# Written by Fang Tian and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# First, define the directory to store the Singularity image
IMG_DIR="$HOME/storage/tmp"

# Then, define the paths
MCR_dir="/apps/matlab/MATLAB_Compiler_Runtime/v95"
input_dir="$CBIG_TESTDATA_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/single_subject"
out_dir=$1 # Your output directory
curr_sub="sub005"
anat=${curr_sub}_FS

# ! Note: the following three paths are paths within the container, which has /input mounted to $input_dir
container_anat_dir="/input/data"
container_fmrinii_dir="/input/docker/sub005_docker.fmrinii"
container_config_file="/input/docker/prepro_docker.config"

curr_dir=$(pwd)
work_dir=$HOME/cluster/

echo $curr_dir
echo $work_dir

if [ ! -e $work_dir ]; then
    mkdir -p $work_dir
fi

cd $work_dir

log_file="${out_dir}/CBIG_preproc_unit_tests_call_fMRI_preproc_singularity.log"

# Define the Singularity command
cmd="source ~/.bashrc && module load singularity/3.11.1 && \
export SINGULARITY_CACHEDIR=$IMG_DIR/.singularity/cache && \
conda activate CBIG_py3 && \
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/CBIG_preproc_run_container.sh \
-c singularity -m $MCR_dir -i $input_dir -o $out_dir \
-sdir $IMG_DIR \
-s $curr_sub -anat_s $anat \
-anat_d $container_anat_dir \
-fmrinii $container_fmrinii_dir \
-config $container_config_file \
-nocleanup | tee -a ${log_file}"

# Submit the job
temp_script_file="${out_dir}/temp_script.sh"
echo '#!/bin/bash' >>${temp_script_file}
echo ${cmd} >>${temp_script_file}
chmod 755 ${temp_script_file}

$CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd "${temp_script_file}" -walltime 2:00:00 -mem 16G \
    -name "CBIG_preproc_unit_tests_call_fMRI_preproc_singularity"
