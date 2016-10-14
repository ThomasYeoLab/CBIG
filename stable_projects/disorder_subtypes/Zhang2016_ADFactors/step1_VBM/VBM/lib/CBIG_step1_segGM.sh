#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -b <brainList> -o <outDir> -t <brainList_tmp> [-q <queue>]
	- brainList		Text file with each line being the path to a brain image; e.g., ~/outputs/VBM/brainList.txt
	- brainList_tmp		(Optional) subset of brainList indicating which brains to include in constructing the study-specific template; if not provided, all will be included
	- queue			(Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs; if not provided, jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":b:o:t:q:" opt; do
	case "${opt}" in
		b) brainList=${OPTARG};;
        	o) outDir=${OPTARG};;
        	t) brainList_tmp=${OPTARG};;
        	q) queue=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${brainList}" ] || [ -z "${outDir}" ] || [ -z "${brainList_tmp}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

echo 'Step 1: GM Segmentation'

mkdir -p ${outDir}
# Will echo outputs into these two files along the way
> ${outDir}GMList.txt
> ${outDir}GMList_tmp.txt
# Relative paths to absolute
outDir=$(readlink ${outDir} -f)/
brainList_tmp=$(readlink ${brainList_tmp} -f)

for brain in `cat ${brainList}`; do
	tmpVar="${brain##*/}" # discard everything before /
	filename="${tmpVar%.nii.gz}" # discard the extension
	if [ -z "${queue}" ]; then
		export outDir brainList_tmp brain filename
		./CBIG_step1_segGM_job.sh
	else
		logDir=${outDir}logs/
		mkdir -p ${logDir}
		qsub -q ${queue} -v outDir=${outDir},brainList_tmp=${brainList_tmp},brain=${brain},filename=${filename} -o ${logDir}${filename}.out -e ${logDir}${filename}.err ./CBIG_step1_segGM_job.sh
	fi
done
./CBIG_waitUntilFinished.sh ${outDir}GMList.txt ${brainList}

echo 'Step 1: GM Segmentation -- Finished'

