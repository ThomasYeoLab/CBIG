#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script projects fsaverage surface data to MNI152 volume space using Freesurfer baseline approaches

###########################################
#Define paths
###########################################

UTILITIES_DIR=$(dirname "$(readlink -f "$0")")/utilities
fs_ver=`cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' | head -c 3`
DEFAULT_MASK=$(dirname "$(dirname "$(readlink -f "$0")")")/bin/liberal_cortex_masks_FS$fs_ver/FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz

###########################################
#Main commands
###########################################
main(){
  #Loop through each input
  input_paths=`cat $input_list`
  for input in $input_paths
  do
    input_name=`basename $input .mat`
    output_prefix=$input_name.fsaverage_to_${template_sub_id}_$fs_type
    #Write input into curvature files
    matlab -nosplash -nojvm -nodesktop -r "load('$input'); \
                                           write_curv('lh.temp.curv', lh_label, 327680); \
                                           write_curv('rh.temp.curv', rh_label, 327680); \
                                           exit"

    #Project input to volume
    cmd="$UTILITIES_DIR/CBIG_RF_${fs_type}_proj_surf2vol.sh $template_sub_dir $template_sub_id temp $output_dir $output_prefix $interp"
    echo $cmd
    eval $cmd

    #Set up file names
    input_lh=$output_dir/lh.$output_prefix.nii.gz
    input_rh=$output_dir/rh.$output_prefix.nii.gz
    output=$output_dir/$output_prefix.nii.gz
    output_seg=$output_dir/seg.$output_prefix.nii.gz

    #Grow the output volume and make seg file
    matlab -nodesktop -nojvm -nosplash -r "addpath('$UTILITIES_DIR/'); \
                                           [output_vol, output_seg_vol] = CBIG_RF_propagate_in_vol('$input_lh', '$input_rh', '$mask'); \
                                           MRIwrite(output_vol, '$output'); \
                                           MRIwrite(output_seg_vol, '$output_seg'); \
                                           exit"

    #Remove intermediate files
    rm lh.temp.curv rh.temp.curv
  done
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -l <input_list> -f <fs_type> -s <template_sub_id> -d <template_sub_dir> -i <interp> -c <cortical_mask> -o <output_dir>

This script projecs input(s) from fsaverage surface to a volumetric atlas space using Affine or MNIsurf/Colinsurf approach.

REQUIRED ARGUMENTS:
	-l <input_list> 	absolute path to file containing input data file names. Each line in the file should contain one file name. Note that the inputs should be in .mat format, containing two variables 'lh_label' and 'rh_label', corresponding to the input data in left and right hemispheres respectively.

OPTIONAL ARGUMENTS:
	-f <fs_type>            type of FS baseline approach to use (Affine or MNIsurf/Colinsurf). For Colinsurf, also call MNIsurf.
				[ default: MNIsurf ]
	-s <template_sub_id>	subject ID of the volumetric template used in recon-all
				[ default: FSL_MNI152_FS4.5.0 ]
	-d <template_sub_dir>	SUBJECTS_DIR of the volumetric template's recon-all results
				[ default: $CBIG_CODE_DIR/data/templates/volume/ ]
	-i <interp> 		interpolation (trilinear or nearest)
				[ default: nearest ]
	-c <cortical_mask> 	absolute path to liberal cortical mask. This should be a binary mask that loosely defines the cortex of the volumetric atlas space. The shape of the mask will determine the shape of the output.
				[ default: $DEFAULT_MASK ]
	-o <output_dir> 	absolute path to output directory
				[ default: $(pwd)/results/projected_fsaverage2vol ]
	-h			display help message

OUTPUTS:
	$0 will create 2 files for each input, corresponding to the projected data in the volumetric atlas sapce and the projected data in segmentation form (with left hemisphere values starting from 0 and right hemisphere values starting from 1000). 
	For example: 
		input_name.fsaverage_to_FSL_MNI152_FS4.5_Affine.nii.gz
		seg.input_name.fsaverage_to_FSL_MNI152_FS4.5_Affine.nii.gz

EXAMPLE:
	$0 -l my_data_list.csv
	$0 -l my_data_list.csv -s SPM_Colin27_FS4.5.0 -n 50

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

###########################################
#Parse arguments
###########################################

#Default parameters
fs_type=MNIsurf
template_sub_id=FSL_MNI152_FS4.5.0
template_sub_dir=$CBIG_CODE_DIR/data/templates/volume/
interp=nearest
mask=$DEFAULT_MASK
output_dir=$(pwd)/results/projected_fsaverage2vol/

#Assign arguments
while getopts "l:f:s:d:i:c:o:h" opt; do
  case $opt in
    l) input_list=${OPTARG} ;;
    f) fs_type=${OPTARG} ;;
    s) template_sub_id=${OPTARG} ;;
    d) template_sub_dir=${OPTARG} ;;
    i) interp=${OPTARG} ;;
    c) mask=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

###########################################
#Check parameters
###########################################

if [ -z $input_list ]; then
  echo "Input list not defined."; 1>&2; exit 1
fi

###########################################
#Other set-ups
###########################################

#Make sure output directory is set up
if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
#Implementation
###########################################

main
