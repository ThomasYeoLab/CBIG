#!/bin/sh

#PBS -l walltime=01:00:0
#PBS -l mem=4gb

bet ${img} ${outDir}${name}_brain.nii.gz -f 0.4 -g -0.1 -B
rm ${outDir}${name}_brain_mask.nii.gz
