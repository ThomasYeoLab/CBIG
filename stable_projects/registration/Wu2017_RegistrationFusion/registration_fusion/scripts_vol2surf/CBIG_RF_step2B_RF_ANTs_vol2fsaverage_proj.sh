#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script projects MNI152 index volumes to a specifiec number of GSP subjects' T1 space, using ANTs transforms, and then to fsaverage surface.

###########################################
#Define paths
###########################################

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
UTILITIES_DIR=$(dirname "$SCRIPT_DIR")/utilities
DEFAULT_GSP_SUBLIST=$(dirname "$(dirname "$SCRIPT_DIR")")/bin/GSP_subjectid.csv
ANTSREG_PREP_SCRIPT=$(dirname "$SCRIPT_DIR")/CBIG_RF_ANTsReg_prep.sh

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

      #Set up variables
      input=$input_dir/${template_type}_$xyz.INDEX.nii.gz
      output_prefix=${xyz}Index_RF_ANTs_${template_type}_to_$sub
      intermediate=$output_dir/index_T1/$output_prefix.nii.gz
      output=$output_dir/index_fsaverage/rh.${output_prefix}_to_fsaverage.nii.gz

      if [ ! -e $output ]; then
        #Create temporarry .nii file for the subject's norm.mgz
        reference=$output_dir/temp_norm_$sub.$xyz.nii
        mri_convert $ind_sub_dir/$sub/mri/norm.mgz $reference

        #Project the index files
        cmd="CBIG_antsApplyReg_vol2vol.sh -i $input -r $reference -d $warp_dir -w $warp_prefix -o $output_dir/index_T1 -p $output_prefix -s inverse -t linear -a $ANTs_dir; rm $reference; \
             $UTILITIES_DIR/CBIG_RF_proj_T1_to_fsaverage.sh $intermediate $output_dir/index_fsaverage $output_prefix $ind_sub_dir ${sub}"

         #Submit a job to PBS scheduler if specified. Otherwise the command is executed directly
         if [ ! -z $queue ]; then
           $UTILITIES_DIR/CBIG_RF_imgRegProj_pbsubmit.sh $queue $output_dir/job_files RF_ANTs_proj 5 1 8 "$cmd"
           sleep $interval
         else
           echo $cmd
           eval $cmd
        fi
      else
        echo "$sub has already been projected to fsaverage."
      fi
    done
  done
}

###########################################
#Function usage & set-up
###########################################

#usage
usage() { echo "
Usage: $0 -p <template_type> -n <num_of_sub> -a <ants_dir> -i <input_dir> -w <warp_dir> -l <ind_sub_list> -g <ind_sub_dir> -o output_dir -q <queue> -t <interval>

This script projects existing x/y/z/ index files in a volumetric atlas space to fsaverage surface as step 2 in RF-ANTs approach. The index files are registered to the volumetric atlas space using ANTs registration results (see $ANTSREG_PREP_SCRIPT about how to run ANTs registration). They are then projected to fsaverage surface.

REQRUIED ARGUMENTS:
	-p <template_type>      type of volumetric template used in step 1, index files creation. See $SCRIPT_DIR/CBIG_RF_step1_make_xyzIndex_volTemplate.sh for more details.

OPTIONAL ARGUMENTS:
	-n <num_of_sub>		number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.
				[ default: 0 ]
	-a <ants_dir> 		absolute path to directory where ANTs is installed
				[ default: $ANTSPATH ]
	-i <input_dir> 		absolute path to input directory. The inputs are the index files created in step 1, i.e. the input directory should be the same as output directory in step 1.
				[ default: $(pwd)/results/index_MNI152 ]
	-w <warp_dir> 		absolute path to ANTs registration results directory. If the resutls were generated using $ANTSREG_PREP_SCRIPT, this directory should be the same as output directory used in the preparation step.
				[ default: $(pwd)/results/antsReg ]
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
	1) index_T1 folder: 3 files will be generated for each subject, corresponding to the x/y/z index files projected to the subject's T1 space. 
	For example: 
		xIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS.nii.gz
                yIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS.nii.gz
                zIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS.nii.gz
	2) index_fsaverage folder: 3 files will be generated for each subject, corresponding to the x/y/z index files projected to fsaverage through that subject. 
	For example: 
		lh.xIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
                rh.xIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
                lh.yIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
                rh.yIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
                lh.zIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
                rh.zIndex_RF_ANTs_FSL_MNI152_FS4.5_to_Sub0001_Ses1_FS_to_fsaverage.nii.gz
	3) job_files folder: contains output and error log files for PBS scheduler (only created if -q has been set)

EXAMPLE:
	$0 -p 'MNI152_orig' -q circ-spool
	$0 -p 'my_template' -i $(pwd)/results/my_index/ -n 50

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
input_dir=$(pwd)/results/index_MNI152/
warp_dir=$(pwd)/results/antsReg/
ind_sub_list=$DEFAULT_GSP_SUBLIST
ind_sub_dir=/mnt/yeogrp/data/GSP_release/
output_dir=$(pwd)/results
interval=10s

#Assign arguments
while getopts "p:n:a:i:w:l:g:q:t:o:h" opt; do
  case $opt in
    p) template_type=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    a) ANTs_dir=${OPTARG} ;;
    i) input_dir=${OPTARG} ;;
    w) warp_dir=${OPTARG} ;;
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

if [ -z $template_type ]; then
  echo "Template type not defined."; 1>&2; exit 1
fi

###########################################
#Other set-ups
###########################################

#Make sure output directory is set up
if [ ! -d "$output_dir/index_T1" ]; then
  echo "Output directory for index in T1 does not exist. Making directory now..."
  mkdir -p $output_dir/index_T1
fi
if [ ! -d "$output_dir/index_fsaverage" ]; then
  echo "Output directory for index in fsaverage does not exist. Making directory now..."
  mkdir -p $output_dir/index_fsaverage
fi

###########################################
#Implementation
###########################################

main


