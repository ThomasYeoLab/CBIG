#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script linearly registers a volume file from its original template's space to its norm.mgz space

###########################################
#Main commands
###########################################
main(){

  base_input=`basename $input .nii.gz`
  #Register the header information
  cmd="mri_vol2vol --mov $template_orig --targ $template_norm --regheader --reg-final $output_dir/$base_input.orig2norm.dat --out $output_dir/$base_input.orig2norm.nii.gz"
  echo $cmd
  eval $cmd

  #Convert warp to FSL format
  cmd="tkregister2 --mov $template_orig --targ $template_norm --reg $output_dir/$base_input.orig2norm.dat --fslregout $output_dir/$base_input.orig2norm.mat --noedit"
  echo $cmd
  eval $cmd

  #Apply registration to input mask
  cmd="flirt -in $input -ref $template_norm -applyxfm -init $output_dir/$base_input.orig2norm.mat -interp nearestneighbour -o $output_dir/$output_prefix.nii.gz -v"
  echo $cmd
  eval $cmd

  #Removei intermediate file
  rm $output_dir/$base_input.orig2norm.*
}

##################################################################
#Function usage
##################################################################

#usage
usage() { echo "
Usage: $0 input template_orig template_norm output_dir output_prefix

This script linearly registers a volume file from its original template's space to its norm.mgz space

	<input>		absolute path to input volume
	<template_orig>	absolute path to the original template
	<template_norm> absolute path to the template's norm.mgz file
	<output_dir>	absolute path to output directory
	<output_prefix>	prefix for output warps. 

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
#Assign input variables
##################################################################
input=$1
template_orig=$2
template_norm=$3
output_dir=$4
output_prefix=$5

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
