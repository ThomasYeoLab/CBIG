#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -p <paramfile> -t <test_fold> -o <outdir>
    This function computes the prediction regression for each cross-validation test folds
    - paramfile     Full path of the .mat file that stores prediction parameters
    - test_fold     Index of the test fold to run the inner-loop cross-validation
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

LF=$outdir/logs/innerloop_${test_fold}.log
mkdir -p $outdir/job_err_out
matlab -nodesktop -nosplash -nodisplay -r " load $paramfile; CBIG_KRR_innerloop_cv_allparams_LITE( $test_fold,\
    param.sub_fold, param.num_inner_folds,param.outdir, param.outstem, param.with_bias, 'none', param.ker_param,\
    param.lambda_set,param.threshold_set, param.metric ); exit; " >> $LF 2>&1