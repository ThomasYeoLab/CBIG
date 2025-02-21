#!/bin/bash
# This script runs CBIG_fMRI_Preproc2016 single subject unit test within a Docker container.
# ! Warning: Since Docker is not installed on the compute nodes, this unit test runs directly on the compiler.
# !          Therefore, please choose a time when the compiler is not busy to run this unit test.
#
# Written by Fang Tian and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

DOCKER_IMG_NAME="thomasyeolab/cbig_fmri_preproc2016:latest"

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

log_file="${out_dir}/CBIG_preproc_unit_tests_call_fMRI_preproc_docker.log"

# Define the Docker command
cmd="source ~/.bashrc && conda activate CBIG_py3 && \
$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/container/CBIG_preproc_run_container.sh \
-c docker -m $MCR_dir -i $input_dir -o $out_dir \
-s $curr_sub -anat_s $anat \
-anat_d $container_anat_dir \
-fmrinii $container_fmrinii_dir \
-config $container_config_file \
-nocleanup | tee -a ${log_file}"

# Execute the command
eval $cmd
