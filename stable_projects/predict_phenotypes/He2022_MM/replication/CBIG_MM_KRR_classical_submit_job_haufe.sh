#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

# check whether have already run
if [[ "$8" == *HCP_diff_roi_pfc.mat ]]; then
    prefix='HCP'
else
    prefix='ukbb'
fi

RESULT_DIR="$2/output_phe_$1/${prefix}_$1_k_$5_rng_num_$9"
FILE=$RESULT_DIR/final_result.mat
if [ -f "$FILE" ]; then
    echo "$FILE exist"
    exit 0
fi

# sleep 5s to reduce job scheduler load
sleep 5s

# remove old log directory if exists, create new one
LOG_DIR="$2/log/${prefix}_krr_$1/fewshot_k_$5"
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}

# submit job
JOB=`$CBIG_SCHEDULER_DIR/qsub -V <<EOJ
#!/bin/bash
#PBS -l ncpus=1
#PBS -l walltime=6:00:00
#PBS -l mem=10gb
#PBS -m ae
#PBS -e ${LOG_DIR}/stderr.txt
#PBS -o ${LOG_DIR}/stdout.txt
matlab -nosplash -nodisplay -nodesktop -r " \
        addpath('$4'); \
        CBIG_MM_KRR_classical('$CBIG_CODE_DIR', '$7', '$3', '$8', '$2', '$9', '$1', '$5', '1'); \
        exit;"
echo ${LOG_DIR}
EOJ`

echo "JobID = ${JOB} for $1 and $5 shot submitted"
