#!/bin/bash
#####
# Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/MRtrix/ \
#        CBIG_DiffProc_tractography_4_del_tractograms.sh $subj_list $output_dir
#
# This script removes tractograms produced by the CBIG MRtrix pipeline to save storage space. A list of 
# subjects and the output directory is required 
# (e.g. /mnt/isilon/CSC1/Yeolab/Data/HCP/HCP_task_derivatives/tractography/iFOD2/output).
#
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

###############
# set up environment
###############
list=$1
output_dir=$2

###############
# delete tractograms
###############
curr_dir=$( pwd )
# move to tmp directory to prevent accidental deletion of files
mkdir -p tmp
cd tmp

# remove all files except for connectomes folder
cat $list | while read subject; do
    sub=$( echo $subject | tr -d '\r' )
    echo $sub
    if [ -e $output_dir/$sub ]; then
        cd $output_dir/$sub
        echo "Deleting tractogram for subj ID: $sub"
        ls | grep -v "connectomes" | xargs rm -r
        cd $curr_dir
    else
        echo "WARNING: Could not find directory"
    fi
done

# remove tmp directory
cd ..
rm -r tmp