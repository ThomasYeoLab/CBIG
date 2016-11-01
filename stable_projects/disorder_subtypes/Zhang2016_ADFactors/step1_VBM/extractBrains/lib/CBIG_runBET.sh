#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -i <imgList> [-n <nameList>] -p <BETParam> -o <outDir> [-q <queue>]
	- imgList	Text file with each line pointing to a NIfTI volume; e.g., ../replicatePNAS/imgList_BETParam_bl/imgList1.txt
	- nameList	(Optional) text file with the same number of lines as imgList. Each line specifies a shorter (or more meaningful) name of the corresponding image. If not provided, output names will be built upon the original filenames
	- BETParam	Shell script, inside which you can customize the parameters for FSL BET; e.g., ./CBIG_BETParam1.sh
	- outDir	Output directory; e.g., ~/outputs/VBM/brains1/
	- queue		(Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs. If not provided, jobs will run serially
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":i:n:p:o:q:" opt; do
	case "${opt}" in
		i) imgList=${OPTARG};;
		n) nameList=${OPTARG};;
        	p) BETParam=${OPTARG};;
        	o) outDir=${OPTARG};;
        	q) queue=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${imgList}" ] || [ -z "${BETParam}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

mkdir -p ${outDir}

lineNo=0
for img in `cat ${imgList}`; do
	lineNo=$((lineNo+1))

	# Get image name
	if [ -z "${nameList}" ]; then
		# Name list not provided -- use original filenames
		tmpVar="${img##*/}" # discard everything before /
		name="${tmpVar%.nii.gz}"; name="${name%.nii}" # discard the extension
	else
		# Name list provided -- just take the corresponding line
		name=`sed "${lineNo}q;d" ${nameList}`
	fi

	# Run or submit BET job
	if [ -z "${queue}" ]; then
		# No cluster -- run jobs serially
		export img name outDir
		${BETParam}
	else
		# Use cluster -- parallelize jobs
		logDir=${outDir}logs/
		mkdir -p ${logDir}
		outDir=$(readlink ${outDir} -f)/
		qsub -q ${queue} -v img=${img},name=${name},outDir=${outDir} -o ${logDir}${name}.out -e ${logDir}${name}.err ${BETParam}  
	fi

	# Note down the path of img for convinience
	echo ${img} > ${outDir}${name}_origPath.txt
done

