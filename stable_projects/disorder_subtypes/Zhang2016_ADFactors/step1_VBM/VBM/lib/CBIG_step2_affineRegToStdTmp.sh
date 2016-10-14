#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -g <GMList> -o <outDir> [-q <queue>]
	- GMList		Text file with each line being the path to a GM image; e.g., ~/outputs/VBM/GM/GMList_tmp.txt
	- outDir		Output directory; e.g., ~/outputs/VBM/affineRegToStdTmp/
	- queue			(Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs; if not provided, jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":g:o:q:" opt; do
	case "${opt}" in
		g) GMList=${OPTARG};;
        	o) outDir=${OPTARG};;
        	q) queue=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${GMList}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

echo 'Step 2: Affine Registration to FSL Standard Template'

mkdir -p ${outDir}
> ${outDir}GMToStdTmpList.txt
# Relative paths to absolute
outDir=$(readlink ${outDir} -f)/

for GM in `cat ${GMList}`; do
	tmpVar="${GM##*/}" # discard everything before /
	filename="${tmpVar%_pve_1.nii.gz}" # discard the suffix
	if [ -z "${queue}" ]; then
		export outDir filename GM
		./CBIG_step2_affineRegToStdTmp_job.sh
	else
		logDir=${outDir}logs/
		mkdir -p ${logDir}
		qsub -q ${queue} -v GM=${GM},filename=${filename},outDir=${outDir} -o ${logDir}${filename}.out -e ${logDir}${filename}.err ./CBIG_step2_affineRegToStdTmp_job.sh
	fi
done
./CBIG_waitUntilFinished.sh ${outDir}GMToStdTmpList.txt ${GMList}

echo 'Step 2: Affine Registration to FSL Standard Template -- Finished'

