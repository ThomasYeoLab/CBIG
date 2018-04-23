#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This function projects x, y, z index volumes of a single subject (after being warped into each subject's T1 space) onto the Freesurfer fsaverage surface space, using FreeSurfer's mri_vol2surf

###########################################
#Main commands
###########################################
main(){
  #Print FreeSurfer version
  echo "FreeSurfer version: `cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' `"

  #Project input from the subject's T1 space onto Freesurfer fsaverage surface space
  export SUBJECTS_DIR=$sub_dir
  for hemi in lh rh
  do
    output=$output_dir/$hemi.${output_prefix}_to_fsaverage.nii.gz
    cmd="mri_vol2surf --mov $input --hemi $hemi --projfrac 0.5 --trgsubject fsaverage --regheader $sub_id --o $output --reshape --interp trilinear"
    echo $cmd
    eval $cmd
  done
}

##################################################################
#Function usage
##################################################################

#usage
usage() { echo "
Usage: $0 input output_dir output_prefix sub_dir sub_id

This script projects an input volume to fsaverage surface using FreeSurfer's mri_vol2surf

	<input>		absolute path to input volume
	<output_dir>	absolute path to output directory
	<output_prefix>	prefix for output warps. For example, the left hemisphere output is named lh.output_prefix_to_fsaverage.nii.gz
	<sub_dir> 	SUBJECTS_DIR of the input's recon-all results
	<sub_id> 	Subject ID of the input used in recon-all
	
" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
#Assign input variables
##################################################################
input=$1
output_dir=$2
output_prefix=$3
sub_dir=$4
sub_id=$5

##################################################################
#Set up output directory
##################################################################

if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
#Implementation
###########################################

main



