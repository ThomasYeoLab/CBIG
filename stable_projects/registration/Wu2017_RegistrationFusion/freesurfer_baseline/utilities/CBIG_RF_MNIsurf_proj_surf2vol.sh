#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This function script projects a single surface in fsaverage space to a volume space using FS-recon approach. It is required to finish running the recon-all results before running this function.

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
    input=$hemi.$input_prefix.curv
    output=$output_dir/$hemi.${output_prefix}.nii.gz
    template=$sub_dir/$sub_id/mri/norm.mgz
  
    #Project to subject's surface space
    cmd="mri_surf2surf --srcsubject fsaverage --srcsurfval $input --hemi $hemi --trgsubject $sub_id --trgsurfval $output --reshape --mapmethod nnf"
    echo $cmd
    eval $cmd

    #Project to subject's T1 space
    cmd="mri_surf2vol --surfval $output --hemi $hemi --subject $sub_id --identity $sub_id --fillribbon --o $output --template $template"
    echo $cmd
    eval $cmd
  done
}

##################################################################
#Function usage
##################################################################

#usage
usage() { echo "
Usage: $0 sub_dir sub_id input_prefix output_dir output_prefix

This script projects an input surface to a target subject's T1 space, using MNIsurf/Colinsurf approach.

	<sub_dir>	SUBJECTS_DIR of the target subject's recon-all results
	<sub_id> 	Subject ID of the target subject used in recon-all
	<input_prefix> 	prefix of input surface. For example, input in left hemisphere would be lh.input_prefix.curv
	<output_dir> 	absolute path to output directory
	<output_prefix> prefix for output warps. For example, the output from left hemispehre will be named lh.output_prefix.nii.gz

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
input_prefix=$3
output_dir=$4
output_prefix=$5

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
#Make sure the specified subject directory has a fsaverage folder
##################################################################
if [ ! -d "$sub_dir/fsaverage" ]; then
  echo "The subject directory does not contain a fsaverage folder. Making a link from default directory..."
  ln -s $SUBJECTS_DIR/fsaverage $sub_dir
fi

###########################################
#Implementation
###########################################

main






