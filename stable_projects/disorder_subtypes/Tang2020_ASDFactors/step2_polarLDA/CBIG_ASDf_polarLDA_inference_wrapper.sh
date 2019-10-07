#!/bin/sh

# Wrapper script to infer factor compositions of new participants with polarLDA model

# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

###################
# Input variables
###################
corpusDir=$1 # document corpus
modelDir=$2 # learned model
outputDir=$3 # your output directory
outputName=$4 # your output name
codeDir=$5 # your code directory
infSettings=$6 # inference parameters setting file
inputDir=$7    # directory of the estimated E(RSFC patterns|Factor) and Pr(Factor|Participant)
factorNum=$8 # number of factors

####################################
# Get run number of final estimate
####################################
ref_run=$(basename $(readlink -f "${inputDir}/k${factorNum}/r*"))
ref_runNum=${ref_run#*r}

#################
# Run inference
#################
csh ${codeDir}/CBIG_ASDf_polarLDA_inf.csh \
-corpus ${corpusDir} \
-model_dir ${modelDir} \
-factor_num ${factorNum} \
-run_num ${ref_runNum} \
-output_dir ${outputDir} \
-output_name ${outputName} \
-infSettings ${infSettings} \
-code_dir ${codeDir}

########################################
# Write normalized factor compositions
########################################
gamma_file=${outputDir}/k${factorNum}r${ref_runNum}_${outputName}-gamma.dat
factorComp_fileName=${outputDir}/${outputName}_k${factorNum}r${ref_runNum}_factorComp.txt
cd ${codeDir}
matlab -nodisplay -nosplash -nodesktop -r \
"clear all;close all;clc; \
CBIG_ASDf_gamma2table('${gamma_file}','${factorComp_fileName}');exit;"


