#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -p <paramfile> -n <N_test_folds> -o <N_test_folds>
    This function computes the prediction regression for each cross-validation test folds
    - paramfile     Full path of the .mat file that stores prediction parameters
    - N_test_folds  Number of cross-validation test folds
    - outdir        Output directory
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":p:n:o:" opt; do
    case "${opt}" in
        p) paramfile=${OPTARG};;
        n) N_test_folds=${OPTARG};;
        o) outdir=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${paramfile}" ] || [ -z "${N_test_folds}" ] || [ -z "${outdir}" ] ; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################
 
LF=$outdir/logs/testloop.log
mkdir -p $outdir/job_err_out
matlab -nodesktop -nosplash -nodisplay -r " load $paramfile; for i = 1:${N_test_folds};\
CBIG_KRR_test_cv_allparams_LITE( i, param.sub_fold, param.outdir, param.outstem, param.with_bias,\
    'none',param.ker_param, param.lambda_set, param.threshold_set );end; exit; " >> $LF 2>&1