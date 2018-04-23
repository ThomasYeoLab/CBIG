#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script runs ANTs registration from subject's T1 space to MNI152 template for a specified number of GSP subjects

###########################################
#Define paths
###########################################

DEFAULT_GSP_SUBLIST=$(dirname "$(dirname "$(readlink -f "$0")")")/bin/GSP_subjectid.csv
UTILITIES_DIR=$(dirname "$(readlink -f "$0")")/utilities

###########################################
#Main commands
###########################################
main(){
  #Get subject names
  if [ $num_sub -gt 0 ]; then
    sub_names=`head -$num_sub $ind_sub_list`
    else
      sub_names=`cat $ind_sub_list`
  fi

  #Loop through each subject
  for sub in $sub_names
  do
    #Create temporarry .nii file for the subject's norm.mgz
    input=$output_dir/antsReg/temp_norm_$sub.nii
    mri_convert $ind_sub_dir/$sub/mri/norm.mgz $input

    #Set up output names
    output_prefix=${sub}_moving_${template_type}_fixed_ants
    output=$output_dir/antsReg/${output_prefix}1InverseWarp.nii.gz

    #Run registration if output does not exist yet
    if [ ! -e $output ]; then
      cmd="CBIG_antsReg_vol2vol.sh -r $template -i $input -d $output_dir/antsReg -p $output_prefix -a $ANTs_dir; rm $input"

      #Submit a job to PBS scheduler if specified. Otherwise the command is executed directly
      if [ ! -z $queue ]; then
        $UTILITIES_DIR/CBIG_RF_imgRegProj_pbsubmit.sh $queue $output_dir/job_files ants_reg 20 1 15 "$cmd"
        sleep $interval
      else
        echo $cmd
        eval $cmd
      fi
    else
      echo "${sub} has already been warped using ANTs"
      rm $input
    fi

  done
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -n <num_of_sub> -a <ants_dir> -m <mni_template> -l <gsp_sub_list> -g <gsp_sub_dir> -q <queue> -t <interval> -o <output_dir>

This script runs ANTs registration between a specified number of individual subjects' T1 space and a volumetric atlas space. The output from this script are used in step 2 in RF-ANTs approach.

REQUIRED ARGUMENTS:
	-p <template_type> 	type of volumetric template to use. Input to this option is also used as prefix of output files.
				Possible options are:
				'MNI152_orig': use FSL_MNI152 1mm template
				'Colin27_orig': use SPM_Colin27 1mm template
				others: use user-defined volumetric template file. In this case, input to this option can be any string; the user is expected to provide a template using "-t" option. Note that volumetric templates themselves are preferred to their normalised volumes, since ANTs registration can take very long to run on normalised volumes (with larger volume dimensions).

OPTIONAL ARGUMENTS:
	-n <num_of_sub>		number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.
				[ default: 0 ]
	-a <ants_dir> 		directory where ANTs is installed
				[ default: $CBIG_ANTS_DIR ]
	-t <template> 		absolute path to user-defined volumetric template file
				[ default: unset ] 	
	-l <ind_sub_list> 	absolute path to a file containing individual subject IDs. Each line in the file should contain one subject ID.
				[ default: $DEFAULT_GSP_SUBLIST ]
	-g <ind_sub_dir> 	SUBJECTS_DIR of individual subjects' recon-all results
				[ default: /mnt/yeogrp/data/GSP_release/ ]
        -o <output_dir>         absolute path to output directory
				[ default: $(pwd)/results ]
	-q <queue> 		for PBS scheduler users, this is equivalent to the -q option for qsub. For example, setting "-q circ-spool" will make the script submit jobs to job scheduler using "qsub -q circ-spool"
				[ default: unset ]
	-i <interval> 		time interval between job submits. For example, the default setting means after each job is submitted, the script 'sleep' for 10 minutes before submitting the next one.
				[ default: 10m ]
	-h			display help message

OUTPUTS:
	$0 will create 2 folders.
	1) antsReg folder: 3 files will be created for each subject, corresponding to the affine transform, forward and backward nonlinear warp from the subject's T1 space to the volumetric atlas space. 
	For example: 
		Sub0001_Ses1_FS_moving_MNI152_orig_fixed_ants0GenericAffine.mat
		Sub0001_Ses1_FS_moving_MNI152_orig_fixed_ants1Warp.nii.gz
		Sub0001_Ses1_FS_moving_MNI152_orig_fixed_ants1InverseWarp.nii.gz
	2) job_files folder: contains output and error log files for PBS scheduler (only created if -q has been set)

EXAMPLE:
	$0 -p 'MNI152_orig' -q circ-spool
	$0 -p 'my_template' -t path/to/my/template.nii.gz -n 50

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

###########################################
#Parse arguments
###########################################

#Default parameters
num_sub=0
ANTs_dir=$CBIG_ANTS_DIR
ind_sub_list=$DEFAULT_GSP_SUBLIST
ind_sub_dir=/mnt/yeogrp/data/GSP_release/
interval=10m
output_dir=$(pwd)/results

#Assign arguments
while getopts "p:n:a:t:l:g:q:i:o:h" opt; do
  case $opt in
    p) template_type=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    a) ANTs_dir=${OPTARG} ;;
    t) template=${OPTARG} ;;
    l) ind_sub_list=${OPTARG} ;;
    g) ind_sub_dir=${OPTARG} ;;
    q) queue=${OPTARG} ;;
    i) interval=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

#Set up default type templates
case $template_type in
  MNI152_orig)
    template=$FSL_DIR/data/standard/MNI152_T1_1mm_brain.nii.gz ;;
  Colin27_orig)
    template=$output_dir/temp_$template_type.nii
    mri_convert $CBIG_CODE_DIR/data/templates/volume/SPM_Colin27_FS4.5.0/mri/orig/001.mgz $template ;;
esac

###########################################
#Check parameters
###########################################

if [ -z $template_type ]; then
  echo "Template type not defined."; 1>&2; exit 1
fi

if [ -z $template ]; then
  echo "User-defined template is not provided."; 1>&2; exit 1
fi

##################################################################
#Set up output directory
##################################################################

if [ ! -d "$output_dir/antsReg" ]; then
  echo "Output directory for registration warp does not exist. Making directory now..."
  mkdir -p $output_dir/antsReg
fi


###########################################
#Implementation
###########################################

main
