#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <input_dir> -f <feature_file> -s <sub_fold> -n <num_score> -o <outdir> -m <model> [-c <cluster>]
    This function computes the predictive-feature matrices for a regression model
    - input_dir     Directory where the regression results are stored
    - feature_file  Full path of the regression feature file
    - sub_fold      Full path of the sub_fold_file that stores the cross-validation split indices
    - num_score     Number of behavioral scores to compute the PFM
    - outdir        Output directory
    - model         regression model. Choose from singleKRR, multiKRR, and LRR
    - cluster       (Optional) if you do not have a cluster, put it as "none", then the computation will
		    run serially (potentially very slow!) for each behavior. If you have a cluster, put
		    your cluster name and the function submit parallel jobs to your cluster
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:f:s:n:o:m:c:" opt; do
    case "${opt}" in
        i) input_dir=${OPTARG};;
        f) feature_file=${OPTARG};;
        s) sub_fold=${OPTARG};;
        n) num_score=${OPTARG};;
        o) outdir=${OPTARG};;
        m) model=${OPTARG};;
        c) cluster=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${input_dir}" ] || [ -z "${feature_file}" ] || [ -z "${sub_fold}" ] || \
     [ -z "${num_score}" ] || [ -z "${model}" ] || [ -z "${outdir}" ] ; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

mkdir -p $outdir/job_err_out
mkdir -p $outdir/logs
LF="$outdir/logs/compute_${model}_PFM.log"
date >> $LF

########## compute the PFM

scripts_dir=`dirname "$(readlink -f "$0")"`

for score_ind in `seq 1 $num_score`
do
    echo $score_ind

    if [ "$cluster" == "none" ];then

        ${scripts_dir}/CBIG_TRBPC_PFM_job.sh $input_dir $feature_file $sub_fold $score_ind $outdir $model

    else

        cmd="${scripts_dir}/CBIG_TRBPC_PFM_job.sh $input_dir $feature_file $sub_fold $score_ind $outdir $model"
        errfile=${outdir}/job_err_out/${model}_PFM_score${score_ind}.err
	outfile=${outdir}/job_err_out/${model}_PFM_score${score_ind}.out
	${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 20:00:0 -mem 24gb -joberr ${errfile} -jobout ${outfile}\
		 -cmd "${cmd}" -name TRBPC_PFM

    fi
done
