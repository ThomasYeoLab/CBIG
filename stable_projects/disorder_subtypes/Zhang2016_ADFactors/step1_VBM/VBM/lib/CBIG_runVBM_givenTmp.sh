#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -b <brainList> -t <tmp> -s <sigma> -o <outDir> [-q <queue>]
	- brainList	Text file with each line being the path to a brain image; e.g., ~/outputs/VBM/brainList.txt
	- tmp		Study-specific GM template
	- sigma		Sigma (in mm) of the gaussian kernel for smoothing
	- outDir	Output directory; e.g., ~/outputs/VBM/
	- queue		(Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs; if not provided, jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":b:t:s:o:q:" opt; do
	case "${opt}" in
		b) brainList=${OPTARG};;
        	t) tmp=${OPTARG};;
        	s) sigma=${OPTARG};;
        	o) outDir=${OPTARG};;
        	q) queue=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${brainList}" ] || [ -z "${tmp}" ] || [ -z "${sigma}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi

###########################################
# Main
###########################################

###### Step 1: GM Segmentation
outDir_step1=${outDir}GM/
if [ -z "${queue}" ]; then
	./CBIG_step1_segGM.sh -b ${brainList} -o ${outDir_step1} -t ${brainList}
else
	./CBIG_step1_segGM.sh -b ${brainList} -o ${outDir_step1} -t ${brainList} -q ${queue}
fi
rm ${outDir_step1}GMList_tmp.txt # this file is meaningless under this context


###### Step 6: Nonlinear Registration to Nonlinear Template
outDir_step6=${outDir}nonlinRegToNonlinTmp/
# For ALL GM images (not just those for template construction)
if [ -z "${queue}" ]; then
	./CBIG_step6_nonlinRegToNonlinTmp.sh -g ${outDir_step1}GMList.txt -t ${tmp} -o ${outDir_step6}
else
	./CBIG_step6_nonlinRegToNonlinTmp.sh -g ${outDir_step1}GMList.txt -t ${tmp} -o ${outDir_step6} -q ${queue}
fi


###### Step 7: Concatenate Registered GM and Generate Binary GM Mask
outDir_step7=${outDir}concat/
./CBIG_step7_concatAndGenMask.sh -g ${outDir_step6}GMToNonlinTmpList.txt -o ${outDir_step7}
rm ${outDir_step7}GMToNonlinTmp_mod_mean_binThr*.nii.gz # this file is meaningless under this context


###### Step 8: Smooth
outDir_step8=${outDir}smoothing/
./CBIG_step8_smooth.sh -v ${outDir_step7}GMToNonlinTmp_mod_4d.nii.gz -s ${sigma} -o ${outDir_step8}

