#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This function projects x, y, z index files of a single subject (which contains coordinates in fsaverage projected to subject's surface space) to the subject's T1 space, using FreeSurfer's mri_surf2vol

###########################################
#Main commands
###########################################
main(){
  #Print FreeSurfer version
  echo "FreeSurfer version: `cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' `"

  #Project x, y, and z index volumes from the subject's surface space to the subject's T1 space
  export SUBJECTS_DIR=$sub_dir
  for hemi in lh rh
  do
    input=$input_dir/$hemi.$input_prefix.index
    output=$output_dir/$hemi.${input_prefix}_T1.nii.gz
    template=$sub_dir/$sub_id/mri/norm.mgz
    if [ ! -e $output ]; then
      cmd="mri_surf2vol --surfval $input --hemi $hemi --subject $sub_id --fillribbon --identity $sub_id --o $output --template $template"
      echo $cmd
      eval $cmd
    else
      echo "The input for $sub_id in $hemi has already been projected to T1 space."
    fi
  done
}

##################################################################
#Function usage
##################################################################

#usage
usage() { echo "
Usage: $0s input_dir input_prefix output_dir sub_dir sub_id

This script projects an input surface to a subject's T1 space, using FreeSurfer's mri_surf2vol

	<input_dir> 	absolute path to input directory
	<input_prefix>	prefix of input surface. For example, the input in left hemisphere should be lh.input_prefix.index
	<output_dir> 	absolute path to output directory
	<sub_dir>	SUBJECTS_DIR of the input's recon-all results
	<sub_id> 	Subject ID of the input used in recon-all

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
#Assign input variables
##################################################################
input_dir=$1
input_prefix=$2
output_dir=$3
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



