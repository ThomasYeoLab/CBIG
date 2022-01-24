#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -b <behav_ind> -s <fold_start> -n <fold_num> -o <outdir>
    This function performs the prediction for cross-validation folds
    - behav_ind     The index of the behavior you want to predict
    - fold_start    A scalar. The function only performs the prediction for cross-validation folds from
    		    folder_start to <fold_start + fold_num - 1>
    - fold_num      A scalar. See fold_start.
    - outdir        Output directory
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":b:s:n:o:" opt; do
    case "${opt}" in
        b) behav_ind=${OPTARG};;
        s) fold_start=${OPTARG};;
        n) fold_num=${OPTARG};;
        o) outdir=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${fold_start}" ] || [ -z "${fold_num}" ] || [ -z "${outdir}" ] || [ -z "${behav_ind}" ]; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

mkdir -p ${outdir}/job_err_out
LF="${outdir}/logs/LRR_behav${behav_ind}_foldStart${fold_start}.txt"
date >> $LF

param_file=$outdir/param.mat
scripts_dir=`dirname "$(readlink -f "$0")"`

matlab -nodesktop -nosplash -nodisplay -r " addpath $scripts_dir; load $param_file; \
CBIG_TRBPC_LRR_fitrlinear_workflow_parallel_folds(param,${behav_ind},${fold_start},${fold_num}); exit; " \
>> ${outdir}/logs/LRR_behav${behav_ind}_foldStart${fold_start}.txt 2>&1