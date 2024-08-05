#!/bin/bash

# Make sure jobs are runnig on gpuserver 4
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
log_dir=$ROOTDIR'/job_logs'
job_prefix="DRB"
cd $ROOTDIR

# submission
checkpoint_path=$ROOTDIR'/checkpoints/unmatch2match/harm_model/DeepResBat/ADNI-AIBL'
cmd="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step5_fit_DeepResBat.sh ADNI-AIBL"
job_name=$job_prefix"_Step5_ADNI-AIBL"

sleep 5s
queue=gpuQ
qsub -V -q ${queue} <<EOJ
#!/bin/sh
### the job queue
#PBS -N $job_name
#PBS -l select=1:ncpus=4:mem=48gb:ngpus=1:host=gpuserver4
#PBS -l walltime=16:00:00
#PBS -e $log_dir"/"$job_name".err"
#PBS -o $log_dir"/"$job_name".out"
    $cmd
EOJ
