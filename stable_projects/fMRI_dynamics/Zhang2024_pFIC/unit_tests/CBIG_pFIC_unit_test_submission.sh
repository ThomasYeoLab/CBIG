#!/bin/bash
# This script is for CBIG-users only. The script submits the pFIC unit-test to the CSCHPC job scheduler.
# It first checks the number of available GPU for each GPU server. For replicability, only gpuserver1/2/3 
# are checked. If it finds a gpuserver with at least 1 available GPU, it will submit the job to that gpuserver.
# Otherwise, it will submit the job to gpuserver1. As I do not fully understand how the job scheduler works
# in details, thus I could not exclude the possibility that the job will be stuck at the waiting stage. If that
# happens, users should consider removing the 'host' option (line 46) and proceed. However,
# note that to pass the unit test, the GPU used must be RTX3090, meaning that only gpuserver1/2/3 should be used.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# IMPORTANT: ssh to the headnode before running this script!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The exit status of the unit-test (PASS or FAIL) can be found by opening ~/cluster/pFIC_unittest.o<job_id>

# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# checking the number of available GPUs of each server
num_available_gpu1=`pbsnodes -a|grep '^gpuserver\|resources_assigned.ngpus' |  \
    awk '/gpuserver1/{getline; if (/resources_assigned.ngpus/) print}' | cut -d'=' -f2`
num_available_gpu2=`pbsnodes -a|grep '^gpuserver\|resources_assigned.ngpus' |  \
    awk '/gpuserver2/{getline; if (/resources_assigned.ngpus/) print}' | cut -d'=' -f2`
num_available_gpu3=`pbsnodes -a|grep '^gpuserver\|resources_assigned.ngpus' |  \
    awk '/gpuserver3/{getline; if (/resources_assigned.ngpus/) print}' | cut -d'=' -f2`

# save the unit-test PASS or FAIL result under this folder
output_msg_dir=/home/`whoami`/cluster/ 

if [ ! -d "$output_msg_dir" ]; then
    mkdir -p "$output_msg_dir"
fi

cd $output_msg_dir
cmd_script="$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test.sh"

if (( num_available_gpu3 < 8 )); then
    node=3
elif (( num_available_gpu2 < 8 )); then
    node=2
elif (( num_available_gpu1 < 8 )); then
    node=1
else
    node=1
fi

echo "Job sent to gpuserver${node}"

qsub -V -q gpuQ -l walltime=1:00:00 -l select=1:mem=5G:ngpus=1:host=gpuserver${node} \
    -N "pFIC_unittest" -m ae -- ${cmd_script}

exit 0
