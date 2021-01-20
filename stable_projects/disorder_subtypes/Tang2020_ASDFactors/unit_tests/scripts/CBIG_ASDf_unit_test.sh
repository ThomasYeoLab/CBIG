#!/bin/sh
#
# Example:
#	sh ./CBIG_ASDf_unit_test.sh ~/storage/Temporaray/CBIG_ASDf_unit_test
#
# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

########################################
# Set paths and create log file
########################################

# code_dir is the directory where Tang2019_ASDFactors codes are located
code_dir=${CBIG_CODE_DIR}/stable_projects/disorder_subtypes/Tang2020_ASDFactors
export MATLABPATH=$CBIG_CODE_DIR/setup
# data_dir is the directory where unit tests data is located
unit_test_dir=$CBIG_TESTDATA_DIR/stable_projects/disorder_subtypes/Tang2020_ASDFactors
data_dir=${unit_test_dir}/data

# ref_dir is the directory where reference results are located
ref_dir=${unit_test_dir}/results

# pass in output_dir
output_dir=$1
mkdir -p ${output_dir}/logs
mkdir -p ${output_dir}/job_err_out

# create log file in output_dir
LF="${output_dir}/logs/CBIG_ASDf_unit_test.log"
touch $LF
echo -e "Start CBIG_ASDf_unit_test\n" >> $LF
date >> $LF
echo -e "\nScript directory = ${code_dir}\n" >> $LF
echo -e "Output directory = ${output_dir}\n" >> $LF
echo -e "Reference directory = ${ref_dir}\n\n" >> $LF

##########################
# Define input variables
##########################
# step1A
corrMat_ASD_est=${data_dir}/corrMat_ASD_est.mat
corrMat_con_est=${data_dir}/corrMat_Con_est.mat
sub_info_file_est=${data_dir}/subInfo_est.csv
behav_score_file_est=${data_dir}/behavior_scores_est.csv
output_dir_step1a=${output_dir}/FC2doc
code_dir_step1a=${code_dir}/step1_FC2doc

# step1B
corrMat_ASD_inf=${data_dir}/corrMat_ASD_inf.mat
corrMat_con_inf=${data_dir}/corrMat_Con_inf.mat
ref_inf=${output_dir_step1a}/step1_output_reg_CN_mean_std.mat
sub_info_file_inf=${data_dir}/subInfo_inf.csv
behav_score_file_inf=${data_dir}/behavior_scores_inf.csv
output_dir_step1b=${output_dir}/FC2doc
code_dir_step1b=${code_dir}/step1_FC2doc

# step2A
corpusDir_est=${output_dir_step1a}/step1_output_dx1.dat # output from previous step
clusterName=circ-spool
code_dir_step2a=${code_dir}/step2_polarLDA
infSettings=${code_dir_step2a}/CBIG_ASDf_polarLDA_infSettings.txt
output_dir_step2a=${output_dir}/estimate
factorNum_est=2
progressFile=${output_dir_step2a}/k${factorNum_est}/progress.txt
numRuns_est=5
run_files=${code_dir_step2a}/run_files_5.txt

# step2B
output_dir_step2b=${output_dir}/visualizeFactors
code_dir_step2b=${code_dir}/step2_polarLDA

# step2C
corpusDir_inf=${output_dir_step1a}/step1_output_inf_dx1.dat
modelDir=${output_dir_step2a}
output_dir_step2c=${output_dir}/inference
outputName_inf="factorComp_inf"
code_dir_step2c=${code_dir}/step2_polarLDA
infSettings=${code_dir_step2c}/CBIG_ASDf_polarLDA_infSettings.txt
factorNum_inf=2

#######################################
# Step1A FC2doc for factor estimation
#######################################

echo -e "Performing step1_FC2doc for factor estimation:\n" >> $LF

cd ${code_dir_step1a}
mkdir -p ${output_dir_step1a}

matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all; \
CBIG_ASDf_FC2doc_estFactors_wrapper('${corrMat_ASD_est}', \
'${corrMat_con_est}','${sub_info_file_est}','${output_dir_step1a}');exit"

echo -e "Step1_FC2doc for factor estimation finished.\n\n" >> $LF

###################################################
# Step1B FC2doc for inferring factor compositions
###################################################

echo -e "Performing step1_FC2doc for inference:\n" >> $LF

cd ${code_dir_step1b}

matlab -nosplash -nodisplay -nodesktop -r \
"clear;clc;close all; \
CBIG_ASDf_FC2doc_infFactorComp_wrapper('${corrMat_ASD_inf}', \
'${corrMat_con_inf}','${ref_inf}','${sub_info_file_inf}','${output_dir_step1b}');exit"

echo -e "step1_FC2doc for inference finished.\n\n" >> $LF

##########################################
# Step2A Estimate factors  with polarLDA
##########################################

echo -e "Performing step2_polarLDA factor estimation:\n" >> $LF

cd ${code_dir_step2a}

mkdir -p ${output_dir_step2a}

######## K = 2 ########
echo -e "--------------K = 2------------\n" >> $LF
echo -e "${numRuns_est} random initializations.\n" >> $LF


${code_dir_step2a}/CBIG_ASDf_polarLDA_est.sh \
    -d ${corpusDir_est} \
    -t ${infSettings} \
    -k ${factorNum_est} \
    -m ${code_dir_step2a} \
    -r ${run_files} \
    -o ${output_dir_step2a} \
    -q CBIG_cluster

## hold until submitted jobs are finished
matlab -nodisplay -nosplash -nodesktop -r \
"clear;clc;close all; \
CBIG_ASDf_checkJobStatus('${progressFile}','${numRuns_est}','600');exit;"

echo -e "K = 2 factor estimation finished.\n\n" >> $LF

#######################################################
# Step2B Get final estimate and visualize the factors
#######################################################

echo -e "Get final estimate and visualize the factors.\n" >> $LF

cd ${code_dir_step2b}
mkdir -p ${output_dir_step2b}

matlab -nodisplay -nosplash -nodesktop -r \
"clear;clc;close all; \
CBIG_ASDf_visualizeFactors('${output_dir_step2a}','${output_dir_step2b}','2');exit;"

echo -e "Plotting factor visualization finished.\n\n" >> $LF

################################################################################
# Step2C Infer factor compositions of ASD participants in the inference sample
################################################################################

echo -e "Performing inference of factor compositions.\n" >> $LF

mkdir -p ${output_dir_step2c}

cd ${code_dir_step2c}

######## K = 2 ########

echo -e "--------------K = 2------------\n" >> $LF

## inference
sh ${code_dir_step2c}/CBIG_ASDf_polarLDA_inference_wrapper.sh \
"${corpusDir_inf}" "${modelDir}" "${output_dir_step2c}" \
"${outputName_inf}" "${code_dir_step2c}" "${infSettings}" \
"${output_dir_step2b}" "${factorNum_inf}"

echo -e "K = 2 inference finished.\n\n" >> $LF

##########################
# End of unit tests
##########################

echo -e "End of CBIG_ASDf_unit_test" >> $LF
date >> $LF
echo -e "================================================\n\n" >> $LF
