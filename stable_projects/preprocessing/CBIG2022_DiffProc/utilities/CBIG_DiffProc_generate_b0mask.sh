#!/bin/sh
#####
# Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_diffusion_processing2022/AMICO/CBIG_DiffProc_runAMICO.sh \
#        --subj_list /path/to/txtfile --dwi_dir /path/to/dwi_images \
#        --output_dir /path/to/output --py_env name_of_AMICO_environment \
#        --mask_output_dir /path/to/stored/b0_brainmask
#
# This function fits the NODDI model for a given list of subjects using the AMICO pipeline.
#
# Written by Leon Ooi.
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

dwi_dir=$1
subj=$2
mask_output_dir=$3

IFS=" " read -a bval_arr < "${dwi_dir}/${subj}.bval" # assumes dwi image and bval has been named in the same way
for idx in ${!bval_arr[@]}; do
    if [[ ${bval_arr[$idx]} == 0 ]]; then
        fslroi ${dwi_dir}/$subj $mask_output_dir/$subj/${subj}_b0 ${idx} 1
        bet $mask_output_dir/$subj/${subj}_b0 $mask_output_dir/$subj/${subj}_bet_b0 -m
        break
    fi;
done