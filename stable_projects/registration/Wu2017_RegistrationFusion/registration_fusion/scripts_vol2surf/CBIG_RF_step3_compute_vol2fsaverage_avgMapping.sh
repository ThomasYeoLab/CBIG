#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script computes average mapping from MNI152 to fsaverage across a specified number of GSP subjects, for a specified RF approach (by supplying the corresponding input prefix)

###########################################
#Set-up
###########################################

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
UTILITIES_DIR=$(dirname "$SCRIPT_DIR")/utilities
BIN_DIR=$(dirname "$(dirname "$(dirname "$(readlink -f "$0")")")")/bin
DEFAULT_GSP_SUBLIST=$(dirname "$(dirname "$SCRIPT_DIR")")/bin/GSP_subjectid.csv
fs_ver=`cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' | head -c 3`

###########################################
#Main commands
###########################################

main(){
  #Create temporary subject list
  temp_sub_list=$output_dir/temp_sub_list.csv
  if [ $num_sub -gt 0 ]; then
    head -$num_sub $ind_sub_list > $temp_sub_list
    else
    cat $ind_sub_list > $temp_sub_list
  fi

  #Set up prefix
  input_prefix=${RF_type}_$template_type
  if [ $num_sub -gt 0 ]; then
    output_prefix=${num_sub}Sub_${input_prefix}_to_fsaverage
    else
    output_prefix=allSub_${input_prefix}_to_fsaverage
  fi

  #Call function to run averaging for the specified ANTs setting
  matlab -nodesktop -nojvm -nosplash -r "addpath('$UTILITIES_DIR'); \
                                         CBIG_RF_compute_vol2surf_avgMapping('$temp_sub_list', '$input_dir', '$input_prefix', '$output_dir', '$output_prefix'); \
                                         exit;"

  #Remove temporary file
  rm $temp_sub_list

  #Unless specified, the average mapping is also copied to the bin
  if [ $copy -ne 0 ]; then
    cp $output_dir/?h.avgMapping_$output_prefix.mat $BIN_DIR/final_warps_FS$fs_ver/
    echo "Warning: Average mapping has been copied to bin. Add -c 0 to command if you do not want to copy it."
  fi
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -p <template_type> -i <input_dir> -r <RF_type> -n <num_of_sub> -l <ind_sub_list> -o <output_dir>

This script computes the average mapping from a volumetric atlas space to fsaverage, across a specified number of individual subject, as step 3 in RF approaches. The averaging is carried out across the x/y/z index files previously projected from the volumetric atlas space to fsaverage through each individual subject.

Note that the average mapping will be copied to final warps folder in the bin by default. Use -c option to specify otherwise.

REQUIRED ARGUMENTS:
	-p <template_type>	type of volumetric template used in index files creation. See $SCRIPT_DIR/CBIG_RF_step1_make_xyzIndex_volTemplate.sh for more details.

OPTIONAL ARGUMENTS:
        -i <input_dir>          absolute path to input directory. The inputs are the index files projected to fsaverage surface through individual subjects in step 2, i.e. the input directory should be the same as output directory in step 2.
				[ default: $(pwd)/results/index_fsaverage ]
        -r <RF_type>            type of RF approaches used in step 2, i.e. RF_M3Z or RF_ANTs
                                [ default: RF_ANTs ]
	-n <num_of_sub>		number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.
				[ default: 0 ]
	-l <ind_sub_list> 	absolute path to a file containing individual subject IDs. Each line in the file should contain one subject ID.
				[ default: $DEFAULT_GSP_SUBLIST ]
	-o <output_dir> 	absolute path to output directory
				[ default: $(pwd)/results/mapping ]
	-c <copy_results>	set this to 0 to prevent the results generated from being copied to the bin folders. By default, the average mapping will be copied to final warps folder in the bin.
				[ default: 1 ]
	-h			display help message

OUTPUTS:
	$0 will create 2 files, corresponding to the average mapping from the volumetric atlas space to left and right hemispheres in fsaverage surface respectively. 
	For example: 
		lh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat
                rh.avgMapping_allSub_RF_ANTs_MNI152_orig_to_fsaverage.mat

EXAMPLE:
	$0 -p 'MNI152_norm'
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
input_dir=$(pwd)/results/index_fsaverage/
num_sub=0
RF_type=RF_ANTs
ind_sub_list=$DEFAULT_GSP_SUBLIST
output_dir=$(pwd)/results/mapping/
copy=1

#Assign arguments
while getopts "r:n:p:l:i:o:c:h" opt; do
  case $opt in
    r) RF_type=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    p) template_type=${OPTARG} ;;
    l) ind_sub_list=${OPTARG} ;;
    i) input_dir=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    c) copy=${OPTARG} ;;
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
if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
#Implementation
###########################################

main




