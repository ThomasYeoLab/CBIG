#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script project xyz index files in a specified number of subjects' surface space to the subject's T1 space, and then to a volumetric atlas space

###########################################
#Define paths
###########################################

UTILITIES_DIR=$(dirname "$(dirname "$(readlink -f "$0")")")/utilities
DEFAULT_GSP_SUBLIST=$(dirname "$(dirname "$(dirname "$(readlink -f "$0")")")")/bin/GSP_subjectid.csv

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
    for xyz in x y z #x/y/z index files
    do

      #Project index to subject's T1 space
      input_prefix=${xyz}Index_fsaverage_to_$sub
      cmd="$UTILITIES_DIR/CBIG_RF_proj_fsaverage_to_T1.sh $input_dir $input_prefix $output_dir/index_T1 $ind_sub_dir ${sub}"
      echo $cmd
      eval $cmd

      #Project index to volumetric atlas space
      for hemi in lh rh
      do
        input=$output_dir/index_T1/$hemi.${xyz}Index_fsaverage_to_${sub}_T1.nii.gz
        output_prefix=$hemi.${xyz}Index_fsaverage_to_${sub}_to_${template_sub_id}_RF_M3Z
        output=$output_dir/index_${template_sub_id}/$output_prefix.nii.gz
        if [ ! -e $output ]; then
          cmd="$UTILITIES_DIR/CBIG_RF_talairachM3zReg_vol2vol.sh $ind_sub_dir ${sub} $input $template_sub_dir $template_sub_id $output_dir/index_${template_sub_id} $output_prefix nearest"

          #Submit a job to PBS scheduler if specified. Otherwise the command is executed directly
          if [ ! -z $queue ]; then
            $UTILITIES_DIR/CBIG_RF_imgRegProj_pbsubmit.sh $queue $output_dir/job_files RF_M3Z_proj 5 1 8 "$cmd"
            sleep $interval
          else
            echo $cmd
            eval $cmd
          fi
        else
          echo "$sub $xyz $hemi has already been projected to $template_sub_id."
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
Usage: $0 -i <input_dir> -s <template_sub_id> -n <num_of_sub> -d <template_sub_dir> -l <ind_sub_list> -g <ind_sub_dir> -q <queue> -t <interval> -o <output_dir>

This script projects existing x/y/z index files in individual subject's surface space to a volumetric atlas space as step 2 in RF-M3Z approach. The index files are projected to each subject's own volumetric space. Then they are registered to the FreeSurfer nonlinear volumetric space and to the volumetric atlas space using FreeSurfer's talairach.m3z transform.

REQUIRED ARGUMENTS:
	-s <template_sub_id> 	subject ID of the volumetric template used in recon-all

OPTIONAL ARGUMENTS:
	-n <num_of_sub>		number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.
				[ default: 0 ]
        -i <input_dir>          absolute path to input directory. The inputs are the index files created in step 1, i.e. the input directory should be the same as output directory in step 1.
				[ default: $(pwd)/results/index_fsaverage ]
	-d <template_sub_dir> 	SUBJECTS_DIR of the volumetric template's recon-all results
				[ default: $CBIG_CODE_DIR/data/templates/volume/ ]
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
	$0 will create 3 folders.
	1) index_T1 folder: 6 files will be generated for each subject, corresponding to the x/y/z index files projected to the subject's T1 space, from left and right hemispheres of the subject's surface respectively. 
	For example: 
		lh.xIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                rh.xIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                lh.yIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                rh.yIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                lh.zIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
                rh.zIndex_fsaverage_to_Sub0001_Ses1_FS_T1.nii.gz
	2) index_\${template_sub_id}: 12 files will be generated for each subject, corresponding to the x/y/z index files from left and right hemispheres, registered to FreeSurfer nonlinear volumetric space and to the volumetric atlas space. 
	For example: 
		lh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.fsNonlinear.nii.gz 
                lh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.nii.gz
                rh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.fsNonlinear.nii.gz
                rh.xIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.nii.gz
                lh.yIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.fsNonlinear.nii.gz
                lh.yIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.nii.gz
                rh.yIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.fsNonlinear.nii.gz
                rh.yIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.nii.gz
                lh.zIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.fsNonlinear.nii.gz
                lh.zIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.nii.gz
                rh.zIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.fsNonlinear.nii.gz
                rh.zIndex_fsaverage_to_Sub0001_Ses1_FS_to_FSL_MNI152_FS4.5.0_RF_M3Z.nii.gz
	3) job_files folder: contains output and error log files for PBS scheduler (only created if -q has been set)

EXAMPLE:
	$0 -s 'FSL_MNI152_FS4.5.0' -q circ-spool
	$0 -s 'SPM_Colin27_FS4.5.0' -n 50

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
input_dir=$(pwd)/results/index_fsaverage/
template_sub_dir=$CBIG_CODE_DIR/data/templates/volume/
ind_sub_list=$DEFAULT_GSP_SUBLIST
ind_sub_dir=/mnt/yeogrp/data/GSP_release/
interval=10s
output_dir=$(pwd)/results

#Assign arguments
while getopts "n:i:s:d:l:g:q:t:o:h" opt; do
  case $opt in
    n) num_sub=${OPTARG} ;;
    i) input_dir=${OPTARG} ;;
    s) template_sub_id=${OPTARG} ;;
    d) template_sub_dir=${OPTARG} ;;
    l) ind_sub_list=${OPTARG} ;;
    g) ind_sub_dir=${OPTARG} ;;
    q) queue=${OPTARG} ;;
    t) interval=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

###########################################
#Check parameters
###########################################

if [ -z $template_sub_id ]; then
  echo "Template subject ID not defined."; 1>&2; exit 1
fi

###########################################
#Implementation
###########################################

main


