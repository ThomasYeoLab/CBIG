#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

# check whether have already run
RESULT_DIR="$2/output_phe_$1/ukbb_$1_rng_num_100"
FILE=$RESULT_DIR/final_result.mat
if [ -f "$FILE" ]; then
    echo "$FILE exist"
    exit 0
fi

# sleep 5s to reduce job scheduler load
sleep 5s

# remove old log directory if exists, create new one
LOG_DIR="$2/log/ukbb_krr_$1"
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}

# submit job
JOB=`$CBIG_SCHEDULER_DIR/qsub -V <<EOJ
#!/bin/bash
#PBS -l ncpus=1
#PBS -l walltime=3:00:00
#PBS -l mem=6gb
#PBS -m ae
#PBS -e ${LOG_DIR}/stderr.txt
#PBS -o ${LOG_DIR}/stdout.txt
matlab -nosplash -nodisplay -nodesktop -r " \
        addpath('$4'); \
        CBIG_MM_KRR_filter('$CBIG_CODE_DIR', '$5', '$3', '$6', '100', '$1', '$2'); \
        exit;"
echo ${LOG_DIR}
EOJ`

echo "JobID = ${JOB} for $1 submitted"
