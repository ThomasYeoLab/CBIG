#!/bin/bash

# Make sure jobs are runnig on gpuserver 1 - 3
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ROOTDIR=$CBIG_CODE_DIR"/stable_projects/harmonization/An2024_DeepResBat"
log_dir=$ROOTDIR'/job_logs'
job_prefix="DRB"
cd $ROOTDIR

# submission
checkpoint_path=$ROOTDIR'/checkpoints/unmatch2match/harm_model/cVAE/'$1
cmd="bash $ROOTDIR/replication/scripts/CBIG_DeepResBat_replica_step3_fit_cVAE.sh "$1
job_name=$job_prefix"_Step3_"$1

sleep 5s
queue=gpuQ
qsub -V -q ${queue} <<EOJ
#!/bin/sh
### the job queue
#PBS -N $job_name
#PBS -l select=1:ncpus=4:mem=48gb:ngpus=1:host=gpuserver1+1:host=gpuserver2+1:host=gpuserver3
#PBS -l walltime=8:00:00
#PBS -e $log_dir"/"$job_name".err"
#PBS -o $log_dir"/"$job_name".out"
    $cmd
EOJ
