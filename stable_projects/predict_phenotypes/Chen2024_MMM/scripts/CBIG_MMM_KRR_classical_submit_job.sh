#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

prefix=$8

RESULT_DIR="$2/output_phe_$1/${prefix}_$1_k_$5_rng_num_100"
FILE=$RESULT_DIR/final_result.mat
if [ -f "$FILE" ]; then
  echo "$FILE exist"
  exit 0
fi

# sleep 5s to reduce job scheduler load
sleep 5s

# remove old log directory if exists, create new one
LOG_DIR="$2/log/${prefix}_krr_$1/fewshot_k_$5"
echo ${LOG_DIR}
echo ${prefix}
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}
echo $CBIG_CODE_DIR, $6, $3, $7, $2, '100', $1, $5, '0'

# submit job
JOB=`$CBIG_SCHEDULER_DIR/qsub -V <<EOJ
#!/bin/bash
#PBS -N MMM_KRR_classical_${prefix}
#PBS -l ncpus=2
#PBS -l walltime=6:00:00
#PBS -l mem=16gb
#PBS -m ae
#PBS -e ${LOG_DIR}/stderr.txt
#PBS -o ${LOG_DIR}/stdout.txt
matlab -nosplash -nodisplay -nodesktop -r " \
        addpath('$4'); \
        CBIG_MMM_KRR_classical('$CBIG_CODE_DIR', '$6', '$3', '$7', '$2', '100', '$1', '$5', '0', '${prefix}'); \
        exit;"
echo ${LOG_DIR}
EOJ`

echo "JobID = ${JOB} for $1 and $5 shot submitted"
