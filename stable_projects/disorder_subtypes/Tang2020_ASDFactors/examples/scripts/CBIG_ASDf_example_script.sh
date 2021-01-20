#!/bin/sh

# This script serves as a toy example of using the code in folder step1_FC2doc 
# and step2_polarLDA in Tang2020_ASDFactors
#
# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#############################################
# Specify output directory & code directory
#############################################

output_dir=$1

#############
# Set paths
#############

# code_dir is the directory where Tang2019_ASDFactors is located
code_dir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Tang2020_ASDFactors

# data_dir is the directory where example data is located
data_dir=${code_dir}/examples/input

# ref_dir is the directory where reference results are located
ref_dir=${code_dir}/examples/results

# create a log file
mkdir -p ${output_dir}/logs
LF=${output_dir}/logs/log_file.txt
touch ${LF}
echo -e "Start CBIG_ASDf_example_script\n" >> ${LF}
date >> ${LF}
echo -e "\nScript directory: ${code_dir}\n" >> ${LF}
echo -e "Input directory : ${data_dir}\n" >> ${LF}
echo -e "Output directory: ${output_dir}\n" >> ${LF}
echo -e "Reference directory: ${ref_dir}\n" >> ${LF}

#######################################
# Step1A FC2doc for factor estimation
#######################################

## step1A input variables
corrMat_ASD_est=${data_dir}/corrMat_ASD_est.mat
corrMat_Con_est=${data_dir}/corrMat_Con_est.mat
sub_info_file_est=${data_dir}/subInfo_est.csv
output_dir_step1a=${output_dir}/FC2doc
code_dir_step1a=${code_dir}/step1_FC2doc

cd ${code_dir_step1a}
mkdir -p ${output_dir_step1a}
echo -e "Performing step1_FC2doc for factor estimation:\n" >> ${LF}

matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all; \
CBIG_ASDf_FC2doc_estFactors_wrapper('${corrMat_ASD_est}', \
'${corrMat_Con_est}', \
'${sub_info_file_est}', \
'${output_dir_step1a}');exit"

echo -e "step1_FC2doc for factor estimation finished.\n\n" >> ${LF}

###################################################
# Step1B FC2doc for inferring factor compositions
###################################################

## step1B input variables
corrMat_ASD_inf=${data_dir}/corrMat_ASD_inf.mat
corrMat_Con_inf=${data_dir}/corrMat_Con_inf.mat
sub_info_file_inf=${data_dir}/subInfo_inf.csv
output_dir_step1b=${output_dir}/FC2doc
code_dir_step1b=${code_dir}/step1_FC2doc
ref_inf=${output_dir_step1a}/step1_output_reg_CN_mean_std.mat

echo -e "Performing step1_FC2doc for inference:\n" >> ${LF}

cd ${code_dir_step1b}

matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all; \
CBIG_ASDf_FC2doc_infFactorComp_wrapper('${corrMat_ASD_inf}', \
'${corrMat_Con_inf}', \
'${ref_inf}', \
'${sub_info_file_inf}', \
'${output_dir_step1b}');exit"

echo -e "step1_FC2doc for inference finished.\n\n" >> ${LF}

##########################################
# Step2A Estimate factors  with polarLDA
##########################################

## step2A input variables
corpusDir_est=${data_dir}/step1_output_dx1.dat #word count document
code_dir_step2a=${code_dir}/step2_polarLDA
infSettings=${code_dir_step2a}/CBIG_ASDf_polarLDA_infSettings.txt
output_dir_step2a=${output_dir}/estimate
clusterName=""
factorNum=2
run_files=${code_dir_step2a}/run_files_1.txt

echo -e "Performing step2_polarLDA factor estimation:\n" >> ${LF}

mkdir -p ${output_dir_step2a}

######## K = 2 ########
echo -e "--------------K = 2------------\n" >> ${LF}
echo -e "Running 1 random initialization.\n" >> ${LF}

cd ${code_dir_step2a}

${code_dir_step2a}/CBIG_ASDf_polarLDA_est.sh \
    -d ${corpusDir_est} \
    -t ${infSettings} \
    -k ${factorNum} \
    -m ${code_dir_step2a} \
    -r ${run_files} \
    -o ${output_dir_step2a}

echo -e "K = 2 factor estimation finished.\n\n" >> ${LF}

#######################################################
# Step2B Get final estimate and visualize the factors
#######################################################

## step2B input variables
output_dir_step2b=${output_dir}/visualizeFactors
code_dir_step2b=${code_dir}/step2_polarLDA

echo -e "Get final estimate and visualize the factors:\n" >> ${LF}

cd ${code_dir_step2b}
mkdir -p ${output_dir_step2b}

matlab -nodisplay -nosplash -nodesktop -r \
"clear;clc;close all;addpath(genpath(fullfile(getenv('CBIG_CODE_DIR'),'utilities','matlab')));\
CBIG_ASDf_visualizeFactors('${output_dir_step2a}','${output_dir_step2b}','2');exit;"

echo -e "Plotting factor visualization finished.\n\n" >> ${LF}

################################################################################
# Step2C Infer factor compositions of ASD participants in the inference sample
################################################################################

## step2C input variables
corpusDir_inf=${data_dir}/step1_output_inf_dx1.dat
modelDir=${output_dir_step2a}
output_dir_step2c=${output_dir}/inference
outputName_inf="factorComp_inf"
code_dir_step2c=${code_dir}/step2_polarLDA
infSettings=${code_dir_step2c}/CBIG_ASDf_polarLDA_infSettings.txt

echo -e "Performing inference of factor compositions.\n" >> ${LF}

mkdir -p ${output_dir_step2c}

cd ${code_dir_step2c}

######## K = 2 ########

echo -e "--------------K = 2------------\n" >> ${LF}

## inference
sh ${code_dir_step2c}/CBIG_ASDf_polarLDA_inference_wrapper.sh \
${corpusDir_inf} \
${modelDir} \
${output_dir_step2c} \
${outputName_inf} \
${code_dir_step2c} \
${infSettings} \
${output_dir_step2b} \
2

echo -e "K = 2 inference finished.\n\n" >> ${LF}

###################################################
# Compare your results to reference results
###################################################

echo -e "Comparing your results to reference results.\n" >> ${LF}

cd ${code_dir}/examples/scripts

matlab -nodisplay -nosplash -nodesktop -r \
"clear; close all;clc; \
CBIG_ASDf_check_example_results('${output_dir}','${LF}');exit;"

echo -e "Checking results finished.\n\n" >> ${LF}

##########################
# End of example codes
##########################

echo -e "End of CBIG_ASDf_example_script\n" >> ${LF}
date >> ${LF}
echo -e "================================================\n\n" >> $LF

