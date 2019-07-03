#!/bin/sh

# Written by Gia H. Ngo and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
K=12

for SEED in `seq 1 100`; do
  ALPHA=100
  ETA=0.01

  JOB=`/apps/sysapps/TORQUE/bin/qsub -l nodes=1:ppn=1,mem=20gb -V -q circ-spool<<EOJ

#!/bin/bash
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=1
  cd $SCRIPT_DIR
  matlab -nosplash -nodisplay -nodesktop << M_PROG
  CBIG_AuthorTopicEM_ReplicateBrainMapResults(${SEED},${K},${ALPHA},${ETA});
M_PROG
EOJ`
done
