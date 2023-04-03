#!/bin/bash
# This script will generate the subject fullpaths for the given subject list
# It checks through all possible runs of each subject if the corresponding fMRI data exists
# A fullpath is then generated and tabulated into a csv file
#
# Input arguments:
# 1) subject_list_path: a text file, each line being a subject index.
# 2) lh_out_file: left hemisphere subject fullpath output.
# 2) rh_out_file: right hemisphere subject fullpath output.
#
# Written by Xiaoxuan Yan and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


subject_list_path=`realpath $1`
lh_out_file=`realpath $2`
rh_out_file=`realpath $3`

declare -a run_arr=("002" "003")

#################
# Generate paths
#################
common_suffix='rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6.nii.gz'

while read -r sub; do
    echo "Now is subject No. $sub"
    lh_cur_line=""
    rh_cur_line=""
    for run in "${run_arr[@]}"; do
        lh_file_dir="${CBIG_GSP_DATA_ROOT_DIR}/${sub}/surf/lh.${sub}_bld${run}_${common_suffix}"
        rh_file_dir="${CBIG_GSP_DATA_ROOT_DIR}/${sub}/surf/rh.${sub}_bld${run}_${common_suffix}"
        if [[ -f "${lh_file_dir}" && -f "${rh_file_dir}" ]]; then
            lh_cur_line="${lh_cur_line} ${lh_file_dir}"
            rh_cur_line="${rh_cur_line} ${rh_file_dir}"
        fi
    done

    # write lh data list
    if [[ -z ${lh_cur_line} ]]; then 
    # if cur subject has no fMRI data, need to write a white space, else matlab readtable would skip this line
        echo "Current subject ${sub} has no lh data."
    else
        echo "${lh_cur_line}" | xargs >> "${lh_out_file}"
    fi

    # write rh data list
    if [[ -z ${rh_cur_line} ]]; then 
    # if cur subject has no fMRI data, need to write a white space, else matlab readtable would skip this line
        echo "Current subject ${sub} has no rh data."
    else
        echo "${rh_cur_line}" | xargs >> "${rh_out_file}"
    fi

done < "${subject_list_path}"
