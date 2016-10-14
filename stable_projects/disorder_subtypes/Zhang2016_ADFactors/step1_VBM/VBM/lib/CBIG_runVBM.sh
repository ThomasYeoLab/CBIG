#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -b <brainList> [-t <brainList_tmp>] -s <sigma> -o <outDir> [-q <queue>]
	- brainList		Text file with each line being the path to a brain image; e.g., ~/outputs/VBM/brainList.txt
	- brainList_tmp		(Optional) subset of brainList indicating which brains to include in constructing the study-specific template; if not provided, all will be included
	- sigma			Sigma (in mm) of the gaussian kernel for smoothing
	- outDir		Output directory; e.g., ~/outputs/VBM/
	- queue			(Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs; if not provided, jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":b:t:s:o:q:" opt; do
	case "${opt}" in
		b) brainList=${OPTARG};;
        	t) brainList_tmp=${OPTARG};;
        	s) sigma=${OPTARG};;
        	o) outDir=${OPTARG};;
        	q) queue=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${brainList}" ] || [ -z "${sigma}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi
if [ -z "${brainList_tmp}" ]; then
	brainList_tmp=${brainList} # use all brains for constructing the template 
fi

###########################################
# Main
###########################################

###### Step 1: GM Segmentation
outDir_step1=${outDir}GM/
if [ -z "${queue}" ]; then
	./CBIG_step1_segGM.sh -b ${brainList} -o ${outDir_step1} -t ${brainList_tmp}
else
	./CBIG_step1_segGM.sh -b ${brainList} -o ${outDir_step1} -t ${brainList_tmp} -q ${queue}
fi


###### Step 2: Affine Registration to FSL Standard Template
outDir_step2=${outDir}affineRegToStdTmp/
if [ -z "${queue}" ]; then
	./CBIG_step2_affineRegToStdTmp.sh -g ${outDir_step1}GMList_tmp.txt -o ${outDir_step2}
else
	./CBIG_step2_affineRegToStdTmp.sh -g ${outDir_step1}GMList_tmp.txt -o ${outDir_step2} -q ${queue}
fi


###### Step 3: Create Study-Specific Affine Template
outDir_step3=${outDir}affineTmp/
./CBIG_step3_createAffineTmp.sh -g ${outDir_step2}GMToStdTmpList.txt -o ${outDir_step3}


###### Step 4: Nonlinear Registration to Affine Template
outDir_step4=${outDir}nonlinRegToAffineTmp/
if [ -z "${queue}" ]; then
	./CBIG_step4_nonlinRegToAffineTmp.sh -g ${outDir_step1}GMList_tmp.txt -t ${outDir_step3}affineTmp.nii.gz -o ${outDir_step4}
else
	./CBIG_step4_nonlinRegToAffineTmp.sh -g ${outDir_step1}GMList_tmp.txt -t ${outDir_step3}affineTmp.nii.gz -o ${outDir_step4} -q ${queue}
fi


###### Step 5: Create Study-Specific Nonlinear Template
outDir_step5=${outDir}nonlinTmp/
./CBIG_step5_createNonlinTmp.sh -g ${outDir_step4}GMToAffineTmpList.txt -o ${outDir_step5}


###### Step 6: Nonlinear Registration to Nonlinear Template
outDir_step6=${outDir}nonlinRegToNonlinTmp/
# For ALL GM images (not just those for template construction)
if [ -z "${queue}" ]; then
	./CBIG_step6_nonlinRegToNonlinTmp.sh -g ${outDir_step1}GMList.txt -t ${outDir_step5}nonlinTmp.nii.gz -o ${outDir_step6}
else
	./CBIG_step6_nonlinRegToNonlinTmp.sh -g ${outDir_step1}GMList.txt -t ${outDir_step5}nonlinTmp.nii.gz -o ${outDir_step6} -q ${queue}
fi


###### Step 7: Concatenate Registered GM and Generate Binary GM Mask
outDir_step7=${outDir}concatAndGenMask/
./CBIG_step7_concatAndGenMask.sh -g ${outDir_step6}GMToNonlinTmpList.txt -o ${outDir_step7}


###### Step 8: Smooth
outDir_step8=${outDir}smoothing/
./CBIG_step8_smooth.sh -v ${outDir_step7}GMToNonlinTmp_mod_4d.nii.gz -s ${sigma} -o ${outDir_step8}

