#!/bin/bash

# Written by Xiuming Zhang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###########################################
# Usage and Reading in Parameters
###########################################

# Usage
usage() { echo "
Usage: $0 -b <brainList> -n <nuisanceVars> -k <K> -o <outDir> [-q <queue>]
	- brainList		Text file with each line being the path to a brain image; e.g., ~/VBM/brainList.txt
	- nuisanceVars		Text file with each line specifying the nuisance variable(s) for the corresponding line in branList
	- K			Number of factors; possible values: 2, 3 and 4, because only these three models are released
	- outDir		Output directory; e.g., ~/outputs/
	- queue			(Optional) if you have a cluster, use it to specify the queue to which you want to qsub these jobs; if not provided, jobs will run serially (potentially very slow!)
" 1>&2; exit 1; }

# Reading in parameters
while getopts ":b:n:k:o:q:" opt; do
	case "${opt}" in
		b) brainList=${OPTARG};;
        	n) nuisanceVars=${OPTARG};;
        	k) K=${OPTARG};;
        	o) outDir=${OPTARG};;
        	q) queue=${OPTARG};;
       		*) usage;;
    	esac
done
shift $((OPTIND-1))
if [ -z "${brainList}" ] || [ -z "${nuisanceVars}" ] || [ -z "${K}" ] || [ -z "${outDir}" ]; then
	echo Missing Parameters!
	usage
fi
if [ "${K}" -ne 2 ] && [ "${K}" -ne 3 ] && [ "${K}" -ne 4 ]; then
	echo K can only be 2, 3 or 4!
	usage
fi

###########################################
# Main
###########################################

###### Step 1: VBM Given GM Template
outDir_step1=${outDir}VBM/
GMTmp=$(readlink ./files_VBM/nonlinTmp.nii.gz -f) # need absolute path
sigma=4.246
cd ../../step1_VBM/VBM/lib # cd there to use code
if [ -z "${queue}" ]; then
	./CBIG_runVBM_givenTmp.sh -b ${brainList} -t ${GMTmp} -s ${sigma} -o ${outDir_step1}
else
	./CBIG_runVBM_givenTmp.sh -b ${brainList} -t ${GMTmp} -s ${sigma} -o ${outDir_step1} -q ${queue}
fi
cd ../../../ADNIDataRelease/inferNew # cd back


###### Step 2: Converting VBM Images to Documents
outDir_step2=${outDir}LDA/
vol4D=${outDir_step1}smoothing/GMToNonlinTmp_mod_4d_s${sigma}.nii.gz
concatOrder=${outDir_step1}concat/GMToNonlinTmp_mod_4d_concatOrder.txt
GMMask=$(readlink ./files_VBM/GMToNonlinTmp_mod_mean_binThr0.05.nii.gz -f) # need absolute path
params=$(readlink ./files_LDA/refParams.mat -f) # need absolute path
# Line order of nuisanceVars corresponds to that of brainList, but not to that of concatOrder, which depends on the completion order of cluster jobs. Hence, we should reorder lines of nuisanceVars so that it and concatOrder correspond
matlab -nodisplay -nosplash -r "brainList = '${brainList}'; nuisanceVars = '${nuisanceVars}'; concatOrder = '${concatOrder}'; CBIG_reorderNuisance(brainList, nuisanceVars, concatOrder); exit;"
# Now the real business
matlab -nodisplay -nosplash -r "cd ../../step2_LDA/lib/; vol4D = '${vol4D}'; concatOrder = '${concatOrder}'; GMMask = '${GMMask}'; params = '${params}'; nuisanceVars = '${nuisanceVars}'; nuisanceVars = [nuisanceVars(1:end-4) '_matchConcatOrder.csv']; outDir = '${outDir_step2}'; CBIG_brain2doc(vol4D, concatOrder, GMMask, params, nuisanceVars, outDir); exit;"


###### Step 3: LDA Inference
docs=${outDir_step2}docs.dat
model=$(readlink ./files_LDA/model_K${K}/final -f) # need absolute path
outName=${outDir_step2}new
cd ../../step2_LDA/lib/ # cd there to use code
./CBIG_LDA_inf.sh -d ${docs} -m ${model} -o ${outName}
cd ../../ADNIDataRelease/inferNew # cd back
