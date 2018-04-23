#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This function script projects a single volume to fsaverage space using MNIsurf approach. It is required to finish running the recon-all results before running this function.

###########################################
#Main commands
###########################################
main(){ 
  #Print FreeSurfer version
  echo "FreeSurfer version: `cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' `"

  #Project the volume to fsaverage surface space
  export SUBJECTS_DIR=$sub_dir
  for hemi in lh rh
  do
    output=$output_dir/$hemi.$output_prefix.nii.gz
    cmd="mri_vol2surf --mov $input --hemi $hemi --projfrac 0.5 --regheader $sub_id --trgsubject fsaverage --o $output --reshape --interp $interp"
    echo $cmd
    eval $cmd
  done
}

##################################################################
#Function usage
##################################################################

#usage
usage() { echo "
Usage: $0 sub_dir sub_id input output_dir output_prefix interp

This script projects an input volume to fsaverage using MNIsurf/Colin27 approach.

	<sub_dir>	SUBJECTS_DIR of the input subject's recon-all results
	<sub_id> 	Subject ID of the input subject used in recon-all
	<input> 	absolute path to input volume
	<output_dir> 	absolute path to output directory
	<output_prefix> prefix for output warps. For example, the output in left hemispehre will be named lh.output_prefix.nii.gz
	<interp>	interpolation (nearest or trilinear)

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
#Set up parameters
##################################################################
sub_dir=$1
sub_id=$2
input=$3
output_dir=$4
output_prefix=$5
interp=$6

##################################################################
#Check if recon-all has been done yet.
##################################################################
if [ ! -d "$sub_dir/$sub_id" ]; then
  echo "Recon-all has not been done yet. Please run recon-all of the input volume before calling this function."
  exit
fi

##################################################################
#Make sure output directory is set up
##################################################################
if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

##################################################################
#Make sure the specified subject directory a fsaverage folder
##################################################################
if [ ! -d "$sub_dir/fsaverage" ]; then
  echo "The subject directory does not contain a fsaverage folder. Making a link from default directory..."
  ln -s $SUBJECTS_DIR/fsaverage $sub_dir
  exit
fi

###########################################
#Implementation
###########################################

main







