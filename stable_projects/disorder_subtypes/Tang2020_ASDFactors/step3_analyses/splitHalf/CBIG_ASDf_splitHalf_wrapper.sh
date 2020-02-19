#!/bin/sh
#
# Estimate factors in two random splits independently
#
# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

output_dir=$1

data_dir=${CBIG_REPDATA_DIR}/stable_projects/disorder_subtypes/Tang2020_ASDFactors/results/splitHalf
corpus_set1=${data_dir}/set1_dx1.dat
corpus_set2=${data_dir}/set1_dx2.dat
K=3

code_dir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Tang2020_ASDFactors/step2_polarLDA

## step2A input variables
output_dir_step2a=${output_dir}/estimate
polarLDA_dir=${code_dir}/polarLDA
infSettings=${code_dir}/CBIG_ASDf_polarLDA_infSettings.txt
progressFile_set1=${output_dir_step2a}/set1/k${K}/progress.txt
progressFile_set2=${output_dir_step2a}/set2/k${K}/progress.txt
run_files=${code_dir}/run_files_50.txt

## step2B input variables
output_dir_step2b=${output_dir}/visualizeFactors


mkdir -p ${output_dir}/job_logs
LF=${output_dir}/job_logs/CBIG_ASDf_step2_polarLDA.log
touch $LF

##############################
# Submit job to circ cluster
##############################

$CBIG_SCHEDULER_DIR/qsub -V -q circ-spool << EOJ

#!/bin/bash
#PBS -S /bin/bash
#PBS -N 'CBIG_ASDf_splitHalf'
#PBS -l walltime=65:00:00
#PBS -l mem=8gb
#PBS -e ${output_dir}/job_logs/CBIG_ASDf_splitHalf.err
#PBS -o ${output_dir}/job_logs/CBIG_ASDf_splitHalf.out

## step2A estimate
${code_dir}/CBIG_ASDf_polarLDA_est.sh \
    -d ${corpus_set1} \
    -t ${infSettings} \
    -k 3 \
    -m ${code_dir} \
    -r ${run_files} \
    -o ${output_dir_step2a}/set1 \
    -q circ-spool

${code_dir}/CBIG_ASDf_polarLDA_est.sh \
    -d ${corpus_set2} \
    -t ${infSettings} \
    -k 3 \
    -m ${code_dir} \
    -r ${run_files} \
    -o ${output_dir_step2a}/set2 \
    -q circ-spool

## hold until submitted jobs are finished
matlab -nodisplay -nosplash -nodesktop -r \
"clear;clc;close all;CBIG_ASDf_checkJobStatus('${progressFile_set1}','50','600');exit;"

matlab -nodisplay -nosplash -nodesktop -r \
"clear;clc;close all;CBIG_ASDf_checkJobStatus('${progressFile_set2}','50','600');exit;"

## step2B visualize factors
cd ${code_dir}
mkdir -p ${output_dir_step2b}

matlab -nodisplay -nosplash -nodesktop -r \
"clear;clc;close all;CBIG_ASDf_visualizeFactors('${output_dir_step2a}', \
'${output_dir_step2b}','${K}');exit;"

EOJ
