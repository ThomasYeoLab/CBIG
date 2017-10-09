#!/bin/sh
# This script is specific to CIRC HPC cluster.
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

ts=$(date +%Y-%m-%d_%H-%M-%S)

USERNAME=`whoami`

# [PLEASE EDIT THIS]
BASE_DIR="/data/users/${USERNAME}/storage/PROJECT_NAME_${ts}"
# The project folder
# The folder name has a time stamp suffix so that different jobs would output to diferent folders by default.
# If you remove the time stamp, make sure that different jobs do not output to the same file,
# otherwise this will affect the cluster storage performance.

# Create the project folder in case it does not exist.
mkdir -p ${BASE_DIR}

# The standard output of the job will be written to this file as the job run.
# This behaviour is different from that of the stdout.txt file specified below
# which is only written once after the job has finished.
# You can use this file to check on the progress of the job.
# Append " | tee -a ${MANUAL_LOG_FILE}" at the end of each command that you want to log the output.
MANUAL_LOG_FILE="${BASE_DIR}/manual_log.txt"

# [PLEASE EDIT THIS]
ARGUMENT_FOR_JOB_SCRIPT=10
# These are extra variables that will be used in the job script.
# These variables can be paths of input files, input arguments, etc...

# CIRC cluster has 39 nodes and each node has 10 processors.
# For example, if you need 4 nodes and 40 processors, change the resource requirement below to nodes=4:ppn=10
# If you need more than 100GB of RAM, you need to request more than 1 node since each node has about 100GB of RAM.
# A rule of thumb is that the number of nodes should be greater than the RAM requested (in GB) divided by 100GB.
# For example, if you need 250GB of RAM, change the requirement to 'mem=250gb' and 'nodes=3:ppn=1'
# The -V argument in the qsub command is to export all of your command prompt enviroment variables into the job enviroment.
JOB=`/apps/sysapps/TORQUE/bin/qsub -V -q circ-spool<<EOJ
#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -m ae
#PBS -e "${BASE_DIR}/stderr.txt"
#PBS -o "${BASE_DIR}/stdout.txt"
    # [PLEASE EDIT THIS] 
    # Call your Matlab functions here
    matlab -nosplash -nodisplay -nodesktop -r " \
        addpath('/data/users/${USERNAME}/'); \
        disp(randn(${ARGUMENT_FOR_JOB_SCRIPT})); \
        exit;" \
    | tee -a ${MANUAL_LOG_FILE}
    # The tee command copies the standard output of the matlab command to the ${MANUAL_LOG_FILE} file.

    echo ${BASE_DIR}

    # Log the job script variables for checking
    echo ${ARGUMENT_FOR_JOB_SCRIPT} | tee -a ${MANUAL_LOG_FILE}
EOJ`

echo "JobID = ${JOB} submitted"
