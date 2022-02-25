#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

# check whether have already run
RESULT_DIR="$1/ukbb_rng_num_$5"
FILE=$RESULT_DIR/final_result.mat
if [ -f "$FILE" ]; then
    echo "$FILE exist"
    exit 0
fi

# sleep 5s to reduce job scheduler load
sleep 5s

# remove old log directory if exists, create new one
LOG_DIR="$1/log/log_rng_$5"
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${RESULT_DIR}

# submit job
JOB=`$CBIG_SCHEDULER_DIR/qsub -V <<EOJ
#!/bin/bash
#PBS -l ncpus=1
#PBS -l walltime=200:00:00
#PBS -l mem=60gb
#PBS -m ae
#PBS -e ${LOG_DIR}/stderr.txt
#PBS -o ${LOG_DIR}/stdout.txt
matlab -nosplash -nodisplay -nodesktop -r " \
        addpath('$3'); \
        CBIG_MM_KRR_MM_base('$CBIG_CODE_DIR', '$7', '$2', '$8', '$1', '$5', '$6', '$9', '${10}'); \
        exit;"
echo ${LOG_DIR}
EOJ`

echo "JobID = ${JOB} for rng $5 and $6 submitted"
