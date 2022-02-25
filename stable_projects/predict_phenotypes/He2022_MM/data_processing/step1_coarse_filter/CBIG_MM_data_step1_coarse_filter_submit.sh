#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

# remove old log directory if exists, create new one
LOG_DIR="log/"
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}

# submit job
JOB=`$CBIG_SCHEDULER_DIR/qsub -V <<EOJ
#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=200:00:00
#PBS -l mem=256gb
#PBS -m ae
#PBS -e ${LOG_DIR}/stderr.txt
#PBS -o ${LOG_DIR}/stdout.txt
cd $CBIG_CODE_DIR/stable_projects/predict_phenotypes/He2022_MM/data_processing/step1_coarse_filter
source activate CBIG_He2022
python coarse_filter.py
echo ${LOG_DIR}
EOJ`

echo "JobID = ${JOB} submitted"
