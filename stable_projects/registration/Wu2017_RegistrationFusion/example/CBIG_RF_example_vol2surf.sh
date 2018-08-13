#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This function runs RF-ANTs vol2surf mapping creation and projection using 1 subject

RF_DIR=$CBIG_CODE_DIR/stable_projects/registration/Wu2017_RegistrationFusion
RF_VOL2SURF_DIR=$RF_DIR/registration_fusion/scripts_vol2surf
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
  if [ ! -d $SUB_DIR/fsaverage ]; then ln -s $FREESURFER_HOME/subjects/fsaverage $SUB_DIR; fi

  # Step 1: create index files in MNI152
  cmd="$RF_VOL2SURF_DIR/CBIG_RF_step1_make_xyzIndex_volTemplate.sh -p Colin27_orig -o $output_dir/index_Colin27"
  echo $cmd
  eval $cmd

  # Step 2: Project index files to fsaverage through each subject
  cmd="$RF_VOL2SURF_DIR/CBIG_RF_step2B_RF_ANTs_vol2fsaverage_proj.sh -p Colin27_orig -n 1 -i $output_dir/index_Colin27 -w $DATA_DIR/CoRR_HNU/subj01 -o $output_dir -a $ANTs_dir -l $sublist -g $SUB_DIR"
  echo $cmd
  eval $cmd

  # Step 3: Generate average mapping
  cmd="$RF_VOL2SURF_DIR/CBIG_RF_step3_compute_vol2fsaverage_avgMapping.sh -p Colin27_orig -i $output_dir/index_fsaverage -n 1 -o $output_dir/mapping -c 0 -l $sublist"
  echo $cmd
  eval $cmd

  # Project a probabilistic map to fsaverage using the average mapping
  input=$DATA_DIR/prob_map_central_sulc.nii.gz
  echo $input > $output_dir/temp.csv
  cmd="$RF_DIR/bin/scripts_final_proj/CBIG_RF_projectVol2fsaverage_batch.sh -l $output_dir/temp.csv -n 1 -d $output_dir/mapping -o $output_dir/projected_vol2fsaverage -p Colin27_orig"
  echo $cmd
  eval $cmd

  # Remove temporary files and directories
  rm $output_dir/temp.csv
  rm -rf $SUB_DIR/fsaverage
  
}

##################################################################
# Function usage
##################################################################

# Usage
usage() { echo "
Usage: $0 -o output_dir

This script generates an RF-ANTs Colin27-to-fsaverage mapping using a single subject and use it to project a probabilistic map. 
The projected map should be compared to $SCRIPT_DIR/example_results/data/lh.projected_central_sulc.nii.gz and $SCRIPT_DIR/example_results/data/rh.projected_central_sulc.nii.gz. 

REQUIRED ARGUMENTS:
	-o <output_dir> 	absolute path to output directory

OPTIONAL ARGUMENTS:
	-a <ANTs_dir>		directory where ANTs is installed 
				[ default: $CBIG_ANTS_DIR ]
	-h			display help message

OUTPUTS:
	$0 will create 5 folders.
	1) index_Colin27 folder: 3 files will be generated, corresponding to the x/y/z index files in Colin27. 
	The file names will be: 
		Colin27_orig_x.INDEX.nii.gz
		Colin27_orig_y.INDEX.nii.gz
		Colin27_orig_z.INDEX.nii.gz
	2) index_T1 folder: 3 files will be generated, corresponding to the x/y/z index files projected to the subject's T1 space. 
	The file names will be:
		xIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1.nii.gz
		yIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1.nii.gz
		zIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1.nii.gz
	3) index_fsaverage folder: 3 files will be generated, corresponding to the x/y/z index files projected to fsaverage through that subject. 
	The file names will be: 
		lh.xIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1_to_fsaverage.nii.gz
		rh.xIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1_to_fsaverage.nii.gz
		lh.yIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1_to_fsaverage.nii.gz
		rh.yIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1_to_fsaverage.nii.gz
		lh.zIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1_to_fsaverage.nii.gz
		rh.zIndex_RF_ANTs_SPM_Colin27_FS4.5_to_subj01_sess1_to_fsaverage.nii.gz
	4) mapping folder: 2 files will be generated, corresponding to the average mapping from Colin27 space to left and right hemispheres in fsaverage surface respectively. 
	The file names will be: 
		lh.avgMapping_5Sub_RF_ANTs_Colin27_orig_to_fsaverage.mat
		rh.avgMapping_5Sub_RF_ANTs_Colin27_orig_to_fsaverage.mat
	5) projected_vol2fsaverage folder: 2 files will be generated, corresponding to the projected data onto fsaverage left and right hemispheres respectively. 
	The file names will be: 
		lh.prob_map_central_sulc.1Sub_RF_ANTs_Colin27_orig_to_fsaverage.nii.gz
		rh.prob_map_central_sulc.1Sub_RF_ANTs_Colin27_orig_to_fsaverage.nii.gz

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



