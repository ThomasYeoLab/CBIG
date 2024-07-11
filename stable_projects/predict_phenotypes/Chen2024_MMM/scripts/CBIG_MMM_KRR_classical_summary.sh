#!/bin/sh
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script is specific to CBIG HPC cluster.

prefix=$1
base_dir=$CBIG_CODE_DIR"/stable_projects/predict_phenotypes/Chen2024_MMM"
rep_dir=$base_dir"/replication"
log_dir=$rep_dir"/log"

# submit job
JOB=`$CBIG_SCHEDULER_DIR/qsub -V <<EOJ
#!/bin/bash
#PBS -N MMM_KRR_classical_summary_${prefix}
#PBS -l ncpus=2
#PBS -l walltime=1:00:00
#PBS -l mem=16gb
#PBS -m ae
#PBS -e ${log_dir}/KRR_classical_summary_${prefix}_err.txt
#PBS -o ${log_dir}/KRR_classical_summary_${prefix}_out.txt
matlab -nosplash -nodisplay -nodesktop -r " \
        CBIG_CODE_DIR = getenv('CBIG_CODE_DIR'); \
        CBIG_REPDATA = getenv('CBIG_REPDATA_DIR'); \
        base_dir = fullfile(CBIG_CODE_DIR, 'stable_projects', 'predict_phenotypes', 'Chen2024_MMM'); \
        code_dir = fullfile(base_dir, 'KRR_CLASSICAL'); \
        cd(code_dir); \
        data_dir = fullfile(base_dir, 'replication', 'output_KRR_classical_${prefix}'); \
        phe_list = fullfile(CBIG_REPDATA, 'stable_projects', \
        'predict_phenotypes', 'Chen2024_MMM', '${prefix}', '${prefix}_phe_list.txt'); \
        n_rng = 100; \
        CBIG_MMM_KRR_classical_summary(CBIG_CODE_DIR, data_dir, phe_list, n_rng, '${prefix}'); \
        exit;"
echo ${log_dir}
EOJ`