#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -p <paramfile> -t <test_fold> -o <test_fold>
    This function computes the prediction regression for each cross-validation test folds
    - paramfile     Full path of the .mat file that stores prediction parameters
    - test_fold     Index of the current cross-validation test fold to predict
    - outdir        Output directory
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":p:t:o:" opt; do
    case "${opt}" in
        p) paramfile=${OPTARG};;
        t) test_fold=${OPTARG};;
        o) outdir=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${paramfile}" ] || [ -z "${test_fold}" ] || [ -z "${outdir}" ] ; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

kernel_folders=$outdir/kernel_folders.mat
LF=$outdir/logs/testloop_${test_fold}.log
mkdir -p $outdir/job_err_out

matlab -nodesktop -nosplash -nodisplay -r " load $paramfile; load $kernel_folders; CBIG_MultiKRR_FindOptTestAcc(\
param.outdir, param.sub_fold,kernel_folders,param.outstem, param.num_inner_folds,${test_fold},\
param.domain,param.group_kernel, param.with_bias, param.threshold, param.ker_param, \
param.metric); exit; " >> $LF 2>&1