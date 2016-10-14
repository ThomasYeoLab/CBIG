#!/bin/sh

#PBS -l walltime=960:00:0
#PBS -l mem=500mb

# Written by CBIG/Xiuming Zhang under under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Print settings for future recall
date
echo Documents: ${docs}
echo K = ${noTopics}
echo Inference settings are:
cat ${curDir}CBIG_LDA_infSettings.txt

# Avoid the same initializations (seeded by time)
sleep $[ ( $RANDOM % 100 )  + 1 ]s

# Run LDA estimation
alpha=$(echo 1/${noTopics} | bc -l) # heuristics
${curDir}lda-c-dist/lda est ${alpha} ${noTopics} ${curDir}CBIG_LDA_infSettings.txt ${docs} random ${runDir}

# When finished, write a line into the progress file
echo Completed >> ${runDir}../../progress.txt
