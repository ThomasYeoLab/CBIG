#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This function script warps an input volume to FreeSurfer's nonlinear volumetric space, and then to the target volume, using talairach.m3z transform file from recon-all results of the input and target volumes.

###########################################
#Main commands
###########################################
main(){
  #Print FreeSurfer version
  echo "FreeSurfer version: `cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' `"

  #Warp the input volumes to the Freesurfer nonlinear volumetric space
  export SUBJECTS_DIR=$input_sub_dir
  output_fs=$output_dir/$output_prefix.fsNonlinear.nii.gz
  cmd="mri_vol2vol --mov $input --s $input_id --targ $FREESURFER_HOME/average/mni305.cor.mgz --m3z talairach.m3z --o $output_fs --no-save-reg --interp $interp"
  echo $cmd
  eval $cmd

  #Warp to the target's space
  export SUBJECTS_DIR=$target_sub_dir
  output=$output_dir/$output_prefix.nii.gz
  cmd="mri_vol2vol --s $target_id --mov $SUBJECTS_DIR/$target_id/mri/norm.mgz --targ $output_fs --m3z talairach.m3z --no-save-reg --o $output --inv-morph --interp $interp"
  echo $cmd
  eval $cmd
}

##################################################################
#Function usage
##################################################################

#usage
usage() { echo "
Usage: $0 input_sub_dir input_id input target_sub_dir target_id output_dir output_prefix interp

This script maps an input volume to a target's normalised T1 image, using talairach.m3z transform and using the FreeSurfer nonlinear volumetric space as an intermediate space. 

	<input_sub_dir> 	SUBJECTS_DIR of the input's recon-all results
	<input_id> 		Subject ID of the input used in recon-all
	<input>			absolute path to input volume
	<target_sub_dir>	SUBJECTS_DIR of the target's recon-all results
	<target_id> 		Subject ID of the target used in recon-all
	<output_dir>		absolute path to output directory
	<output_prefix>		prefix for output warps. The output is named output_prefix.nii.gz
	<interp>		interpolation (trilin, nearest, cubic)

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
#Set up parameters
##################################################################

input_sub_dir=$1
input_id=$2
input=$3
target_sub_dir=$4
target_id=$5
output_dir=$6
output_prefix=$7
interp=$8

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



