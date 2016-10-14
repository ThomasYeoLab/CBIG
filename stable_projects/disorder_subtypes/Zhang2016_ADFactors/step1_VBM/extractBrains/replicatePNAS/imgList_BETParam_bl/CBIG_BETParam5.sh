#!/bin/sh

#PBS -l walltime=01:00:0
#PBS -l mem=4gb

standard_space_roi ${img} ${outDir}${name}_roi.nii.gz -b
bet ${outDir}${name}_roi.nii.gz ${outDir}${name}_brain.nii.gz -f 0.4
rm ${outDir}${name}_roi.nii.gz # remove intermediate files
