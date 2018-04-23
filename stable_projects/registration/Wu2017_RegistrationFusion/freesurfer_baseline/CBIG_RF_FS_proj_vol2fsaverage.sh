#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script projects MNI152 volume to fsaverage surface space using Freesurfer baseline approaches

###########################################
#Define paths
###########################################

UTILITIES_DIR=$(dirname "$(readlink -f "$0")")/utilities

###########################################
#Main commands
###########################################

main(){
  #Loop through each input
  input_paths=`cat $input_list`
  for input in $input_paths
  do
    input_name=`basename $input .gz`
    input_name=`basename $input_name .nii`
    output_prefix=$input_name.${fs_type}_${template_sub_id}_to_fsaverage
    cmd="$UTILITIES_DIR/CBIG_RF_${fs_type}_proj_vol2surf.sh $template_sub_dir $template_sub_id $input $output_dir $output_prefix $interp"
    echo $cmd
    eval $cmd
  done
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -l <input_list> -f <fs_type> -s <template_sub_id> -d <template_sub_dir> -i <interp> -o <output_dir>

This script projects input(s) from a volumetric atlas space to fsaverage surface using Affine or MNIsurf/Colinsurf approach.

REQUIRED ARGUMENTS:
	-l <input_list> 	absolute path to file containing input data file names. Each line in the file should contain one file name. Note that inputs should be in the format of .nii or .nii.gz

OPTIONAL ARGUMENTS:
        -f <fs_type>            type of FS baseline approach to use (Affine or MNIsurf). For Colinsurf, also call MNIsurf.
				[ default: MNIsurf ]
	-s <template_sub_id>	subject ID of the volumetric template used in recon-all
				[ default: FSL_MNI152_FS4.5.0 ]
	-d <template_sub_dir>	SUBJECTS_DIR of the volumetric template's recon-all results
				[ default: $CBIG_CODE_DIR/data/templates/volume/ ]
	-i <interp> 		interpolation (trilinear or nearest)
				[ default: trilinear ]
	-o <output_dir> 	absolute path to output directory
				[ default: $(pwd)/results/projected_vol2fsaverage ]
	-h			display help message

OUTPUTS:
	$0 will create 2 files for each input, corresponding to the projected data onto fsaverage left and right hemispheres respectively. 
	For example: 
		lh.input_name.Affine_FSL_MNI152_FS4.5.0_to_fsaverage.nii.gz
                rh.input_name.Affine_FSL_MNI152_FS4.5.0_to_fsaverage.nii.gz

EXAMPLE:
	$0 -l my_data_list.csv
	$0 -l my_data_list.csv -f Affine -s 'SPM_Colin27_FS4.5.0'

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
interp=trilinear
output_dir=$(pwd)/results/projected_vol2fsaverage

#Assign arguments
while getopts "f:l:s:d:i:o:h" opt; do
  case $opt in
    f) fs_type=${OPTARG} ;;
    l) input_list=${OPTARG} ;;
    s) template_sub_id=${OPTARG} ;;
    d) template_sub_dir=${OPTARG} ;;
    i) interp=${OPTARG} ;;
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
