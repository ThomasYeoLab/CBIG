#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script projects xyz index files in a specified number of subjects' surface space to the subjects' T1 space, and then to the volumetric atlas space

###########################################
#Define paths
###########################################

UTILITIES_DIR=$(dirname "$(dirname "$(readlink -f "$0")")")/utilities
DEFAULT_GSP_SUBLIST=$(dirname "$(dirname "$(dirname "$(readlink -f "$0")")")")/bin/GSP_subjectid.csv
ANTSREG_PREP_SCRIPT=$(dirname "$(dirname "$(readlink -f "$0")")")/CBIG_RF_ANTsReg_prep.sh

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
    warp_prefix=${sub}_moving_${template_type}_fixed_ants
    for xyz in x y z #x/y/z index files
    do

      #Project index to subject's T1 space
      input_prefix=${xyz}Index_fsaverage_to_$sub
      cmd="$UTILITIES_DIR/CBIG_RF_proj_fsaverage_to_T1.sh $input_dir $input_prefix $output_dir/index_T1 $ind_sub_dir $sub"
      echo $cmd
      eval $cmd

      #Project index to volumetric atlas space
      for hemi in lh rh
      do
        input=$output_dir/index_T1/$hemi.${xyz}Index_fsaverage_to_${sub}_T1.nii.gz     
        output_prefix=$hemi.${xyz}Index_fsaverage_to_${sub}_to_${template_sub_id}_RF_ANTs
        output=$output_dir/index_${template_sub_id}/$output_prefix.nii.gz
        if [ ! -e $output ]; then
          cmd="CBIG_antsApplyReg_vol2vol.sh -i $input -r $template -d $warp_dir -w $warp_prefix -o $output_dir/index_$template_sub_id -p $output_prefix -s forward -t NearestNeighbor -a $ANTs_dir; \
               $UTILITIES_DIR/CBIG_RF_convert_orig2norm.sh $output $template $template_norm $output_dir/index_$template_sub_id $output_prefix"

          #Submit a job to PBS scheduler if specified. Otherwise the command is executed directly
          if [ ! -z $queue ]; then
            $UTILITIES_DIR/CBIG_RF_imgRegProj_pbsubmit.sh $queue $output_dir/job_files RF_ANTs_proj 5 1 8 "$cmd"
            sleep $interval
          else
            echo $cmd
            eval $cmd
          fi
        else
          echo "$sub $xyz index file in $hemi has already been projected to $template_sub_id."
        fi
      done
    done
  done
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -p <template_type> -n <num_of_sub> -a <ants_dir> -i <input_dir> -w <warp_dir> -t <template> -l <ind_sub_list> -g <ind_sub_dir> -o <output_dir> -q <queue> -t <interval>

This script projects existing x/y/z index files in individual subject's surface space to a volumetric atlas space as step 2 in RF-ANTs approach. The index files are projected to each subject's own volumetric space. Then they are registered to the the volumetric atlas space using ANTs registration results (see $ANTSREG_PREP_SCRIPT about how to run ANTs registration).

REQUIRED ARGUMENTS:
	-p <template_type> 	type of volumetric template to use.
				Possible options are:
				'MNI152_orig': use FSL_MNI152 1mm template
				'Colin27_orig': use SPM_Colin27 1mm template
				others: use user-defined volumetric template file. In this case, input to this option an be any string; the user is expected to provide a template using "-t" option.
	-s <template_sub_id> 	subject ID of the volumetric template used in recon-all. Input to this option is used as prefix of output files.

OPTIONAL ARGUMENTS:
	-n <num_of_sub>		number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.
				[ default: 0 ]
	-a <ants_dir> 		absolute path to directory where ANTs is installed
				[ default: $ANTSPATH ]
	-i <input_dir> 		absolute path to input directory. The inputs are the index files created in step 1, i.e. the input directory should be the same as output directory in step 1.
				[ default: $(pwd)/results/index_fsaverage ]
	-w <warp_dir> 		absolute path to ANTs registration results directory. If the resutls were generated using $ANTSREG_PREP_SCRIPT, this directory should be the same as output directory used in the preparation step.
				[ default: $(pwd)/results/antsReg ]
	-t <template> 		absolute path to user-defined volumetric template file
				[ default: unset ]
	-l <ind_sub_list> 	absolute path to a file containing individual subject IDs. Each line in the file should contain one subject ID.
				[ default: $DEFAULT_GSP_SUBLIST ]
	-g <ind_sub_dir> 	SUBJECTS_DIR of individual subjects' recon-all results
				[ default: /mnt/yeogrp/data/GSP_release/ ]
        -o <output_dir>         absolute path to output directory
				[ default: $(pwd)/results/ ]
	-q <queue> 		for PBS scheduler users, this is equivalent to the -q option for qsub. For example, setting "-q circ-spool" will make the script submit jobs to job scheduler using "qsub -q circ-spool"
				[ default: unset ]
	-t <interval> 		time interval between job submits. For example, the default setting means after each job is submitted, the script 'sleep' for 10 seconds before submitting the next one.
				[ default: 10s ]
	-h			display help message

OUTPUTS:
	$0 will create 3 folders:
	1) index_T1 folder: 6 files will be generated for each subject, corresponding to the x/y/z index files projected to the subject's T1 space, from left and right hemispheres of the subject's surface respectively. 
	For example: 
		lh.xIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                rh.xIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                lh.yIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                rh.yIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                lh.zIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                rh.zIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
	2) index_\${template_sub_id}: 6 files will be generated for each subject, corresponding to the x/y/z index files from left and right hemispheres, registered to the volumetric atlas space. 
	For example: 
		lh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_ANTs.nii.gz
                rh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_ANTs.nii.gz
                lh.yIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_ANTs.nii.gz
                rh.yIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_ANTs.nii.gz
                lh.zIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_ANTs.nii.gz
                rh.zIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_ANTs.nii.gz
	3) job_files folder: contains output and error log files for PBS scheduler (only created if -q has been set)

EXAMPLE:
	$0 -p 'MNI152_orig' -s 'FSL_MNI152_FS4.5.0' -q circ-spool
	$0 -t 'my_template' -s 'my_template_FS' -t path/to/my/template.nii.gz -n 50

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
ANTs_dir=$ANTSPATH
input_dir=$(pwd)/results/index_fsaverage/
warp_dir=$(pwd)/results/antsReg/
ind_sub_list=$DEFAULT_GSP_SUBLIST
ind_sub_dir=/mnt/yeogrp/data/GSP_release
output_dir=$(pwd)/results
interval=10s

#Assign arguments
while getopts "s:p:n:a:i:w:t:l:g:q:t:o:h" opt; do
  case $opt in
    s) template_sub_id=${OPTARG} ;;
    p) template_type=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    a) ANTs_dir=${OPTARG} ;;
    i) input_dir=${OPTARG} ;;
    w) warp_dir=${OPTARG} ;;
    t) template=${OPTARG} ;;
    l) ind_sub_list=${OPTARG} ;;
    g) ind_sub_dir=${OPTARG} ;;
    q) queue=${OPTARG} ;;
    t) interval=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

#Set up default type templates
case $template_type in
  MNI152_orig)
    template=$FSL_DIR/data/standard/MNI152_T1_1mm_brain.nii.gz
    template_norm=$CBIG_CODE_DIR/data/templates/volume/FSL_MNI152_FS4.5.0/mri/norm.nii.gz ;;
  Colin27_orig)
    template=$output_dir/temp_$template_type.nii
    mri_convert $CBIG_CODE_DIR/data/templates/volume/SPM_Colin27_FS4.5.0/mri/orig/001.mgz $template
    template_norm=$output_dir/temp_${template_type}_norm.nii
    mri_convert $CBIG_CODE_DIR/data/templates/volume/SPM_Colin27_FS4.5.0/mri/norm.mgz $template_norm;;
esac

###########################################
#Check parameters
###########################################

if [ -z $template_type ]; then
  echo "Template type not defined."; 1>&2; exit 1
fi

if [ -z $template_sub_id ]; then
  echo "Template subject ID not defined."; 1>&2; exit 1
fi

###########################################
#Implementation
###########################################

main


