#!/bin/sh
# CBIG_IndCBM_generate_example_list.sh <output_dir>
# This function generate 6 file lists as input for the example. Each list contains the directories of the fmri files 
# either on the surface or in the volume. 
# Input:
#   output_dir:  Path of output folder. 
# Output:
#   6 list files will be saved under <output_dir>/list as ?h_sub?_sess?.txt
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

output_dir=`realpath ${1}`

echo "Generating example lists..."

list_dir="${output_dir}/list"
if [ ! -d ${list_dir} ];then
    mkdir -p ${list_dir}
fi

sub_num=2
sess_num=1
run_num=1

for ((i=1;i<=$((sub_num));i++))
do
    surf_input_dir="$CBIG_CODE_DIR/data/example_data/CoRR_HNU"
    vol_input_dir="$CBIG_CODE_DIR/stable_projects/brain_parcellation/Xue2021_IndCerebellum/examples/input/vol/sub${i}"
    lh_list_file="${list_dir}/lh_sub${i}_list.txt"
    rh_list_file="${list_dir}/rh_sub${i}_list.txt"
    vol_list_file="${list_dir}/vol_sub${i}_list.txt"

    surf_files=$(find ${surf_input_dir} -name "*subj0${i}*fs5*")
    for filename in $surf_files
    do
        if [[ ${filename} =~ "lh" ]]; then
            echo ${filename} >> ${lh_list_file}
        else
            echo ${filename} >> ${rh_list_file}
        fi
    done

    vol_files=$(ls $vol_input_dir)
    for filename in $vol_files
    do
        echo ${vol_input_dir}/${filename} >> ${vol_list_file}
    done
done

MSHBM_dir="${output_dir}/MSHBM/data_list"
if [ ! -d ${MSHBM_dir}/fMRI_list ];then
    mkdir -p ${MSHBM_dir}/fMRI_list
fi

for sess in {1,2}; do
    for sub in {1,2}; do
        # fMRI data
        lh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/\
lh.subj0${sub}_sess${sess}_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz"

        echo $lh_fmri >> ${MSHBM_dir}/fMRI_list/lh_sub${sub}_sess${sess}.txt

        rh_fmri="$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj0$sub/subj0${sub}_sess${sess}/surf/\
rh.subj0${sub}_sess${sess}_bld002_rest_skip4_stc_mc_residc_interp_FDRMS0.2_DVARS50_bp_0.009_0.08_fs6_sm6_fs5.nii.gz"

        echo $rh_fmri >> ${MSHBM_dir}/fMRI_list/rh_sub${sub}_sess${sess}.txt
    done
done

echo "Finish generating example lists."

