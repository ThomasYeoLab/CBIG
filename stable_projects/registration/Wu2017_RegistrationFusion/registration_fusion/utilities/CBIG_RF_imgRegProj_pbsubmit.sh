#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This functions creates a temporary pbs file in the work_dir to submit, according to the input parametres, so that the outputs of the job reside in the work_dir
#Note that input cmd should be inside double quotes if there are white spaces in the command

##################################################################
#Function usage
##################################################################
if [ $# -ne 7 ]; then
  echo "usage: $0 queue work_dir name wall_hr ppn mem "cmd""
  exit
fi

##################################################################
#Assign variables
##################################################################
queue=$1
work_dir=$2
name=$3
wall_hr=$4
ppn=$5
mem=$6
cmd=$7
curr_dir=`pwd`

##################################################################
#Make sure work_dir exists
##################################################################
if [ ! -e $work_dir ]; then
  mkdir -p $work_dir
fi

##################################################################
#Create and submit the temporary pbs file in work_dir
##################################################################
cd $work_dir
echo "#!/bin/bash" >> pbs_submit.pbs
echo "#PBS -S /bin/bash" >> pbs_submit.pbs
echo "#PBS -N $name" >> pbs_submit.pbs
echo "#PBS -q $queue" >> pbs_submit.pbs
echo "#PBS -l walltime=${wall_hr}:00:00" >> pbs_submit.pbs
echo "#PBS -l nodes=1:ppn=$ppn" >> pbs_submit.pbs
echo "#PBS -l mem=${mem}gb" >> pbs_submit.pbs
echo "cd ${curr_dir}; date; ${cmd}; date" >> pbs_submit.pbs
qsub pbs_submit.pbs

##################################################################
#Delete the temporary file and go back to original directory
##################################################################
rm pbs_submit.pbs
cd $curr_dir
