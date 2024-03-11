#!/bin/bash
# This script is specifically used for submitting a job to the gpuQ from MATLAB
# Written by Shaoshi Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# save the unit-test PASS or FAIL result under this folder
output_msg_dir=/home/`whoami`/cluster/ 

if [ ! -d "$output_msg_dir" ]; then
    mkdir -p "$output_msg_dir"
fi

cd $output_msg_dir
cmd_script="$CBIG_CODE_DIR/stable_projects/fMRI_dynamics/Zhang2024_pFIC/unit_tests/CBIG_pFIC_unit_test_matlab.sh"

ssh headnode "qsub -V -q gpuQ -l walltime=1:00:00 -l select=1:mem=5G:ngpus=1:host=gpuserver3 \
    -N "pFIC_unittest" -m ae -- ${cmd_script}"

exit 0
