#!/bin/sh
# CBIG_IndCBM_create_MSHBM_list.sh <surf_list_dir1> <surf_list_dir2> ... <output_dir>
# This function generate file list for the MSHBM model in replication. The folder structure is the same as MSHBM model.
# Input:
#   surf_list_dir?: Surface file list of subject ?. 
#   output_dir:     Path of output folder. 
# Output:
#   Following the folder structure of MSHBM model, data lists will be saved under <output_dir>/data_list/fMRI_list.   
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

para_num=$#
sub_num=$((para_num-1))
sub_list=($@)
unset sub_list[${sub_num}]
output_dir=${!para_num}

sess_num=100
sid=$1
set=$2

echo "Generating data list for ${sub_num} subjects..."

output_dir="${output_dir}/data_list/fMRI_list"
if [ ! -d ${output_dir} ];then
    mkdir -p ${output_dir}
fi

hemilist="lh rh"
for ((n=0;n<${sub_num};n++))
do
    for hemi in ${hemilist}
    do
        list_dir=${sub_list[$n]}
        list_file=`find ${list_dir} -name ${hemi}*`
        run_list=`cat $list_file`
        i=0
        sess_name=""
        for run_file in ${run_list}
        do
            run_name=`basename ${run_file}`
            run_name=${run_name:3:22}
            if [ "${run_name:0:15}"x != "${sess_name}"x ];then
                i=$((i+1))
                if [ $i -gt ${sess_num} ];then
                    break
                fi
                sess_name=${run_name:0:15}
                outfile="${output_dir}/${hemi}_sub$((n+1))_sess${i}.txt"
                if [ -f ${outfile} ];then
                    rm -f ${outfile}
                fi
            fi
            echo ${run_file} >> ${outfile}
        done
    done
    echo "${i} sessions generated."
done

