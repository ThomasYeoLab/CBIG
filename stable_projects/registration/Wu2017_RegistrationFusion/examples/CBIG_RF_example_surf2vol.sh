#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This function runs RF-ANTs vol2surf mapping creation and projection using 1 subject

RF_DIR=$CBIG_CODE_DIR/stable_projects/registration/Wu2017_RegistrationFusion
RF_SURF2VOL_DIR=$RF_DIR/registration_fusion/scripts_surf2vol
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
DATA_DIR=$CBIG_CODE_DIR/data/example_data/RegistrationFusion_example_data
SUB_DIR=$CBIG_CODE_DIR/data/example_data/CoRR_HNU/subj01_FS

###########################################
# Main commands
###########################################
main(){

  # set up for input subject
  sublist=$output_dir/temp_sublist.csv
  echo "subj01_sess1_FS" > $sublist
  #if [ ! -d $SUB_DIR/fsaverage ]; then ln -s $FREESURFER_HOME/subjects/fsaverage $SUB_DIR; fi

  # Step 1: create fsaverage index files in each subject's surface space
  cmd="$RF_SURF2VOL_DIR/CBIG_RF_step1_make_xyzIndex_fsaverage.sh -n 1 -o $output_dir/index_fsaverage -l $sublist -g $SUB_DIR"
  echo $cmd
  eval $cmd

  # Step 2: Project index files to MNI152
  cmd="$RF_SURF2VOL_DIR/CBIG_RF_step2B_RF_ANTs_fsaverage2vol_proj.sh -p Colin27_orig -s SPM_Colin27_FS4.5.0 -n 1 -i $output_dir/index_fsaverage -w $DATA_DIR/CoRR_HNU/subj01 -o $output_dir -a $ANTs_dir -l $sublist -g $SUB_DIR"
  echo $cmd
  eval $cmd

  # Step 3: Generate average mapping
  cmd="$RF_SURF2VOL_DIR/CBIG_RF_step3_compute_fsaverage2vol_avgMapping.sh -s SPM_Colin27_FS4.5.0 -i $output_dir/index_SPM_Colin27_FS4.5.0 -n 1 -o $output_dir -c 0 -l $sublist"
  echo $cmd
  eval $cmd

  # Project a probabilistic map to fsaverage using the average mapping
  input=$DATA_DIR/surface_parcel.mat
  echo $input > $output_dir/temp.csv
  cmd="$RF_DIR/bin/scripts_final_proj/CBIG_RF_projectfsaverage2Vol_batch.sh -l $output_dir/temp.csv -n 1 -d $output_dir/mapping -o $output_dir/projected_fsaverage2vol -s SPM_Colin27_FS4.5.0 -m $RF_DIR/bin/liberal_cortex_masks_FS5.3/SPM_Colin27_FS4.5.0_cortex_estimate.nii.gz"
  echo $cmd
  eval $cmd

  # Remove temporary files and directories
  rm $output_dir/temp.csv
  
}

##################################################################
# Function usage
##################################################################

# Usage
usage() { echo "
Usage: $0 -o output_dir

This script generates an RF-ANTs fsaverage-to-Colin27 mapping using a single subject and use it to project a surface parcellation. 
The projected map should be compared to example_results/projected_surface_parcel.nii.gz. 

REQUIRED ARGUMENTS:
	-o <output_dir> 	absolute path to output directory

OPTIONAL ARGUMENTS:
	-a <ANTs_dir>		directory where ANTs is installed 
				[ default: $CBIG_ANTS_DIR ]
	-h			display help message

OUTPUTS:
	$0 will create 5 folders.
	1) index_fsaverage folder: 6 files will be generated, corresponding to the x/y/z index files in the subject's surface space. 
	The file names will be:  
		lh.xIndex_fsaverage_to_subj01_sess1.index
		rh.xIndex_fsaverage_to_subj01_sess1.index
		lh.yIndex_fsaverage_to_subj01_sess1.index
		rh.yIndex_fsaverage_to_subj01_sess1.index
		lh.zIndex_fsaverage_to_subj01_sess1.index
		rh.zIndex_fsaverage_to_subj01_sess1.index
	2) index_T1 folder: 6 files will be generated, corresponding to the x/y/z index files projected to the subject's T1 space, from left and right hemispheres of the subject's surface respectively. 
	The file names will be: 
		lh.xIndex_fsaverage_to_subj01_sess1_T1.nii.gz
		rh.xIndex_fsaverage_to_subj01_sess1_T1.nii.gz
		lh.yIndex_fsaverage_to_subj01_sess1_T1.nii.gz
		rh.yIndex_fsaverage_to_subj01_sess1_T1.nii.gz
		lh.zIndex_fsaverage_to_subj01_sess1_T1.nii.gz
		rh.zIndex_fsaverage_to_subj01_sess1_T1.nii.gz
	3) index_SPM_Colin27_FS4.5.0: 6 files will be generated, corresponding to the x/y/z index files from left and right hemispheres, registered to the volumetric atlas space. 
	The file names will be: 
		lh.xIndex_fsaverage_to_subj01_sess1_to_SPM_Colin27_FS4.5.0_RF_ANTs.nii.gz
		rh.xIndex_fsaverage_to_subj01_sess1_to_SPM_Colin27_FS4.5.0_RF_ANTs.nii.gz
		lh.yIndex_fsaverage_to_subj01_sess1_to_SPM_Colin27_FS4.5.0_RF_ANTs.nii.gz
		rh.yIndex_fsaverage_to_subj01_sess1_to_SPM_Colin27_FS4.5.0_RF_ANTs.nii.gz
		lh.zIndex_fsaverage_to_subj01_sess1_to_SPM_Colin27_FS4.5.0_RF_ANTs.nii.gz
		rh.zIndex_fsaverage_to_subj01_sess1_to_SPM_Colin27_FS4.5.0_RF_ANTs.nii.gz
	4) mapping folder: corresponding to the average mapping from the fsaverage surface to Colin27 space and the count map in Colin27 space. 
	The file names will be:
		1Sub_fsaverage_to_SPM_Colin27_FS4.5.0_RF_ANTs_avgMapping.prop.mat
		1Sub_fsaverage_to_SPM_Colin27_FS4.5.0_RF_ANTs.count.mat
	5) projected_fsaverage2vol folder: 2 files will be generated, corresponding to the projected data in the volumetric atlas sapce and the projected data in segmentation form (with left hemisphere values starting from 0 and right hemisphere values starting from 1000). 
	The file names will be:
		prob_map_central_sulc.1Sub_fsaverage_to_SPM_Colin27_FS4.5_RF_ANTs.nii.gz
		seg.prob_map_central_sulc.1Sub_fsaverage_to_SPM_Colin27_FS4.5_RF_ANTs.nii.gz

EXAMPLE:
	$0 -o ~/unit_test_results

" 1>&2; exit 1; }

# Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
# Assign input variables
##################################################################

# Default parameter
ANTs_dir=$CBIG_ANTS_DIR

# Assign parameter
while getopts "o:h" opt; do
  case $opt in
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

##################################################################
# Check parameter
##################################################################

if [ -z $output_dir ]; then
  echo "Output directory not defined."; 1>&2; exit 1
fi

##################################################################
# Disable multi-threading
##################################################################

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

##################################################################
# Set up output directory
##################################################################

if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
# Implementation
###########################################

main



