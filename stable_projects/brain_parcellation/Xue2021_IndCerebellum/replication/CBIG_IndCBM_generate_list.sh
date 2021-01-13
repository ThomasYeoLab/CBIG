#!/bin/sh
# CBIG_IndCBM_generate_list.sh <MNI_input_dir> <surf_input_dir> <output_dir>
# This function generate file list for replication for a given subject. 
# Data will be split into discovery set and replication set by session.
# 3 set of lists will be generated: all, discovery, replication.
# Input:
#   MNI_input_dir:  Volume input for all runs.
#   surf_input_dir: Surface input for all runs. Should include left hemisphere and right hemisphere.
#   output_dir:     Path of output folder. 
# Output:
#   5 lists will be saved under <output_dir>/<subj_name>/<set_name> for each subject, each set. 
#   List files lh_<set_name>_list.txt, rh_<set_name>_list.txt, MNI_<set_name>_list.txt contain the fMRI runs. 
#   runs_<set_name>_list.txt contains the name of each run and sess_<set_name>_list.txt contains the name of each 
#   session.
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

MNI_input_dir=$1
surf_input_dir=$2
output_dir=$3

echo "Generating lists from ${MNI_input_dir} and ${surf_input_dir}"
echo "Split by sessions. Discovery set: odd sessions. Replication set: even sessions."

# Find name of all runs
MNI_files=$(ls $MNI_input_dir)
for filename in $MNI_files
do
    run_name=`basename ${filename}`
    run_name=${run_name:0:22}
    run_list="${run_list} ${run_name}"
done

function create_list_dir(){
    # Remove files if exist
    set_name=$1
    list_dir="${output_dir}/${set_name}"
    if [ ! -d ${list_dir} ];then
        mkdir -p ${list_dir}
    fi
    run_list_file="${list_dir}/runs_${set_name}_list.txt"
    if [ -f ${run_list_file} ];then
        rm -f ${run_list_file}
    fi
    sess_list_file="${list_dir}/sess_${set_name}_list.txt"
    if [ -f ${sess_list_file} ];then
        rm -f ${sess_list_file}
    fi
    lh_list_file="${list_dir}/lh_${set_name}_list.txt"
    if [ -f ${lh_list_file} ];then
        rm -f ${lh_list_file}
    fi
    rh_list_file="${list_dir}/rh_${set_name}_list.txt"
    if [ -f ${rh_list_file} ];then
        rm -f ${rh_list_file}
    fi
    MNI_list_file="${list_dir}/MNI_${set_name}_list.txt"
    if [ -f ${MNI_list_file} ];then
        rm -f ${MNI_list_file}
    fi
}

# Generate lists for all runs
set_name="all"
create_list_dir ${set_name}

sess_name=""
numr=0
nums=0
for run_name in $run_list
do
    numr=$((numr+1))
    echo ${run_name} >> ${run_list_file}
    if [ "${run_name:0:15}"x != "${sess_name}"x ];then
        nums=$((nums+1))
        sess_name=${run_name:0:15}
        echo ${sess_name} >> ${sess_list_file}
    fi
    file=$(find ${surf_input_dir} -name "lh.${run_name}_reorient_skip_mc_unwarp_anat_resid_bpss_fsaverage6_sm2.nii.gz")
    echo ${file} >> ${lh_list_file}
    file=$(find ${surf_input_dir} -name "rh.${run_name}_reorient_skip_mc_unwarp_anat_resid_bpss_fsaverage6_sm2.nii.gz")
    echo ${file} >> ${rh_list_file}
    file=$(find ${MNI_input_dir} -name "${run_name}*")
    echo ${file} >> ${MNI_list_file}
done
echo "${nums} sessions, ${numr} runs in total."
echo "Lists saved at ${list_dir}"
allrun_file=${run_list_file}
allsess_file=${sess_list_file}

# Generate lists for discovery set
set_name="discovery"
create_list_dir ${set_name}

i=0
numr=0
nums=0
allrun_list=`cat $allrun_file`
allsess_list=`cat $allsess_file`
for sess_name in $allsess_list
do
    if [ $(($i%2)) == 0 ];then
        nums=$((nums+1))
        echo ${sess_name} >> ${sess_list_file}
        for run_name in $allrun_list
        do
            if [ "${run_name:0:15}"x == "${sess_name}"x ];then
                numr=$((numr+1))
                echo ${run_name} >> ${run_list_file}
                file=$(find ${surf_input_dir} -name \
                    "lh.${run_name}_reorient_skip_mc_unwarp_anat_resid_bpss_fsaverage6_sm2.nii.gz")
                echo ${file} >> ${lh_list_file}
                file=$(find ${surf_input_dir} -name \
                    "rh.${run_name}_reorient_skip_mc_unwarp_anat_resid_bpss_fsaverage6_sm2.nii.gz")
                echo ${file} >> ${rh_list_file}
                file=$(find ${MNI_input_dir} -name "${run_name}*")
                echo ${file} >> ${MNI_list_file}
            fi
        done
    fi
    i=$((i+1))
done
echo "${nums} sessions, ${numr} runs for discovery set."
echo "Lists saved at ${list_dir}"

# Generate lists for replication set
set_name="replication"
create_list_dir ${set_name}

i=0
numr=0
nums=0
allrun_list=`cat $allrun_file`
allsess_list=`cat $allsess_file`
for sess_name in $allsess_list
do
    if [ $(($i%2)) != 0 ];then
        nums=$((nums+1))
        echo ${sess_name} >> ${sess_list_file}
        for run_name in $allrun_list
        do
            if [ "${run_name:0:15}"x == "${sess_name}"x ];then
                numr=$((numr+1))
                echo ${run_name} >> ${run_list_file}
                file=$(find ${surf_input_dir} -name \
                    "lh.${run_name}_reorient_skip_mc_unwarp_anat_resid_bpss_fsaverage6_sm2.nii.gz")
                echo ${file} >> ${lh_list_file}
                file=$(find ${surf_input_dir} -name \
                    "rh.${run_name}_reorient_skip_mc_unwarp_anat_resid_bpss_fsaverage6_sm2.nii.gz")
                echo ${file} >> ${rh_list_file}
                file=$(find ${MNI_input_dir} -name "${run_name}*")
                echo ${file} >> ${MNI_list_file}
            fi
        done
    fi
    i=$((i+1))
done
echo "${nums} sessions, ${numr} runs for replication set."
echo "Lists saved at ${list_dir}"

