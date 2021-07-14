#!/bin/bash

# Written by Jianzhong Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <input_dir> -b <num_behav> -s <outstem> -n <N_per_job> -t <total_num> -o <outdir> -m <model> 
-g <group> [-c <cluster>]
    This function performs the permutation test workflow for kernel regression
    - input_dir     Directory where the kernel regression results are stored.
    - num_behav     Number of behavior scores in the regression model
    - outstem       A string. The outstem you used when you perform the multi-kernel regression.
    - N_per_job     Number of permutations per paralleled job
    - total_num     Total number of permutations
    - outdir        Output directory
    - model         regression model. Choose from singleKRR and multiKRR
    - group         A text file containing the group ID (e.g. site) of all subjects. Permutation will 
    		    only be done within a group
    - cluster       (Optional) if you do not have a cluster, put it as "none", then the computation will
		    run serially (potentially very slow!) for each behavior. If you have a cluster, put
		    your cluster name and the function submit parallel jobs to your cluster
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:b:s:n:t:o:m:g:c:" opt; do
    case "${opt}" in
        i) input_dir=${OPTARG};;
        b) num_behav=${OPTARG};;
        s) outstem=${OPTARG};;
        n) N_per_job=${OPTARG};;
        t) total_num=${OPTARG};;
        o) outdir=${OPTARG};;
        m) model=${OPTARG};;
        g) group=${OPTARG};;
        c) cluster=${OPTARG};;
        *) usage;;
    esac
done
shift $((OPTIND-1))
if [ -z "${input_dir}" ] || [ -z "${num_behav}" ] || [ -z "${outstem}" ] || [ -z "${model}" ] || \
     [ -z "${N_per_job}" ] || [ -z "${total_num}" ] || [ -z "${outdir}" ] ; then
    echo Missing Parameters!
    usage
fi

###########################################
# Main
###########################################

mkdir -p ${outdir}/job_err_out
mkdir -p ${outdir}/logs
LF="${outdir}/logs/perm_input.log"
date >> $LF

########## perform permutation test

scripts_dir=`dirname "$(readlink -f "$0")"`
jobs=0

for score in `seq 1 ${num_behav}`
do
    for perm_start in `seq 1 ${N_per_job} ${total_num}`
    do
        if [ "$cluster" == "none" ];then
            $scripts_dir/CBIG_TRBPC_KRR_perm_job.sh ${model} ${input_dir} ${outstem} ${score} ${perm_start} \
            ${N_per_job} ${group} ${outdir}
        else
            errfile=${outdir}/job_err_out/${model}_perm_score${score}_perm_start${perm_start}.err
            outfile=${outdir}/job_err_out/${model}_perm_score${score}_perm_start${perm_start}.out
            cmd="$scripts_dir/CBIG_TRBPC_KRR_perm_job.sh"
            cmd="${cmd} ${model} ${input_dir} ${outstem} ${score} ${perm_start} ${N_per_job} ${group} ${outdir}"
            ${CBIG_CODE_DIR}/setup/CBIG_pbsubmit -walltime 12:00:0 -mem 16gb -joberr ${errfile} -jobout ${outfile} \
            -cmd "${cmd}" -name ${model}_PFM
        fi
    done
done
