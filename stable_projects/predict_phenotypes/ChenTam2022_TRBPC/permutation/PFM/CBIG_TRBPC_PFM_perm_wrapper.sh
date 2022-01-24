#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <input_dir> -f <feature_file> -s <outstem> -n <perm_num> -t <total_num> -o <outdir> -g <group> 
          [-c <cluster>]
    This function performs the permutation test workflow for the predictive-feature matrices
    - input_dir     Directory where the multi-kernel regression results are stored.
    - feature_file  Full path of the file storing a a cell array. The cell array contains the
                    full path of all the feature files that are required to calculate the 
                    kernels.
    - outstem       A string. The outstem you used when you perform the multi-kernel regression.
    - perm_num      Number of permutations per paralleled job
    - total_num     Total number of permutations
    - outdir        Output directory
    - group         A text file containing the group ID (e.g. site) of all subjects. Permutation will 
    		    only be done within a group
    - cluster       (Optional) if you do not have a cluster, put it as "none", then the computation will
		    run serially (potentially very slow!) for each behavior. If you have a cluster, put
		    your cluster name and the function submit parallel jobs to your cluster
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:f:s:n:t:o:g:c:" opt; do
    case "${opt}" in
        i) input_dir=${OPTARG};;
        f) feature_file=${OPTARG};;
        s) outstem=${OPTARG};;
        n) perm_num=${OPTARG};;
        t) total_num=${OPTARG};;
        o) outdir=${OPTARG};;
        g) group=${OPTARG};;
        c) cluster=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${input_dir}" ] || [ -z "${feature_file}" ] || [ -z "${outstem}" ] || \
     [ -z "${perm_num}" ] || [ -z "${total_num}" ] || [ -z "${outdir}" ] ; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

mkdir -p ${outdir}/job_err_out
mkdir -p ${outdir}/logs
LF="${outdir}/logs/PFM_perm.log"
date >> $LF

########## prepare inputs for the permutation
scripts_dir=`dirname "$(readlink -f "$0")"`

matlab -nodesktop -nosplash -nodisplay -r " addpath $scripts_dir; CBIG_TRBPC_prepare_PFM_perm_inputs( \
   '$feature_file', '$input_dir', '$outdir'); exit; " >> $LF 2>&1
   
########## run permutation jobs
y_file=$outdir/y_regress.mat
FC_mean_file=$outdir/FC_network_mean.mat
krr_folds=$outdir/folds.mat

jobs=0
for perm_seed_start in `seq 1 $perm_num $total_num`
do
    if [ "$cluster" == "none" ];then
    
        $scripts_dir/CBIG_TRBPC_PFM_perm_job.sh $FC_mean_file $y_file $krr_folds $perm_seed_start \
        $perm_num $outdir $group
    else
        errfile=${outdir}/job_err_out/PFM_perm_start_${perm_seed_start}.err
	outfile=${outdir}/job_err_out/PFM_perm_start_${perm_seed_start}.out
	cmd="$scripts_dir/CBIG_TRBPC_PFM_perm_job.sh"
	cmd="$cmd  $FC_mean_file $y_file $krr_folds $perm_seed_start $perm_num $outdir $group"
	${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 12:00:0 -mem 16gb -joberr ${errfile} -jobout ${outfile} \
            -cmd "${cmd}" -name PFM_perm
    fi
done
