#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This function runs the CBIG_make_xyzIndex_surf function to make index volumes for each GSP subject (x/y/z coordinates of vertex in fsaverage registered to subject's surface space)

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
    for hemi in lh rh #left and right hemisphere
    do

      #set up directory and file names
      reg_dir=$ind_sub_dir/${sub}/surf
      output_prefix=fsaverage_to_${sub}
      output=$output_dir/$hemi.${xyz}Index_${output_prefix}.index

      #Run index file creation if the output does not exist yet
      if [ ! -e $output ]; then
        matlab -nodesktop -nosplash -nojvm -r "addpath('$UTILITIES_DIR'); CBIG_RF_make_xyzIndex_surf('$reg_dir', '$output_dir', '$output_prefix'); exit"
      else
        echo "Index file for $sub in $hemi has already been created."
      fi
    done
  done
}

###########################################
#Function usage & set-up
###########################################

#usage
usage() { echo "
Usage: $0 -l <ind_sub_list> -g <ind_sub_dir> -n <num_of_sub> -o <output_dir>

This script creates fsaveraeg x/y/z index files as step 1 in RF approaches. Recon-all process should have been carried out for each individual subject before this step, where surface registration between subject and fsaverage is done. In an index file, each vertex in an individual subject is assigned value based on the x/y/z coordinate of its corresponding vertex in fsaverage surface.

REQUIRED ARGUMENTS:
        -n <num_of_sub>         number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.

OPTIONAL ARGUMENTS:
	-l <ind_sub_list> 	absolute path to a file containing individual subject IDs. Each line in the file should contain one subject ID.
				[ default: $DEFAULT_GSP_SUBLIST ]
	-g <ind_sub_dir> 	SUBJECTS_DIR of individual subjects' recon-all results
				[ default: /mnt/yeogrp/data/GSP_release ]
	-o <output_dir>		absolute path to output directory
				[ default: $(pwd)/results/index_fsaverage ]
	-h			display help message

OUTPUTS:
	$0 will create 6 files for each subject, corresponding to the x/y/z index files in the subject's surface space. 
	For example: 
		lh.xIndex_fsaverage_to_Sub0001_Ses1_FS.index
                rh.xIndex_fsaverage_to_Sub0001_Ses1_FS.index
                lh.yIndex_fsaverage_to_Sub0001_Ses1_FS.index
                rh.yIndex_fsaverage_to_Sub0001_Ses1_FS.index
                lh.zIndex_fsaverage_to_Sub0001_Ses1_FS.index
                rh.zIndex_fsaverage_to_Sub0001_Ses1_FS.index

EXAMPLE:
	$0 -n 0
	$0 -n 50 -g my_sub_dir

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

###########################################
#Parse arguments
###########################################

#Default parameters
ind_sub_list=$DEFAULT_GSP_SUBLIST
ind_sub_dir=/mnt/yeogrp/data/GSP_release
output_dir=$(pwd)/results/index_fsaverage

#Assign arguments
while getopts "l:g:n:o:h" opt; do
  case $opt in
    l) ind_sub_list=${OPTARG} ;;
    g) ind_sub_dir=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

###########################################
#Check parameters
###########################################

if [ -z $num_sub ]; then
  echo "Number of subjects not defined."; 1>&2; exit 1
fi

###########################################
#Other set-ups
###########################################

#Make sure output directory is set up
if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
#Implementation
###########################################

main


