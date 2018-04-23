#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script computes average mapping from fsaverage to the volumetric atlas space across a specified number of GSP subjects for the specified approach approach (by supplying the corresponding input prefix)

###########################################
#Set-up
###########################################

UTILITIES_DIR=$(dirname "$(dirname "$(readlink -f "$0")")")/utilities
BIN_DIR=$(dirname "$(dirname "$(dirname "$(readlink -f "$0")")")")/bin
DEFAULT_GSP_SUBLIST=$BIN_DIR/GSP_subjectid.csv
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
  input_prefix=${template_sub_id}_$RF_type
  if [ $num_sub -gt 0 ]; then
    output_prefix=${num_sub}Sub_fsaverage_to_$input_prefix
    else
    output_prefix=allSub_fsaverage_to_$input_prefix
    num_sub=`cat $temp_sub_list | wc -l`
  fi  

  #Call function to compute average mapping
  matlab -nodesktop -nosplash -nojvm -r "addpath('$UTILITIES_DIR'); \
                                         CBIG_RF_compute_surf2vol_avgMapping('$temp_sub_list', '$input_dir', '$input_prefix', '$output_dir/mapping', '$output_prefix'); \
                                         exit"

  #Call function to create cortical mask
  matlab -nodesktop -nosplash -nojvm -r "addpath('$UTILITIES_DIR'); \
                                         CBIG_RF_make_cortexMask('$template_sub_id', '$output_dir/mapping/${output_prefix}_count.mat', $num_sub, '$output_dir/mask'); \
                                         exit"

  #Remove temporary file
  rm $temp_sub_list

  #Unless specified, the average mapping is also copied to the bin
  if [ $copy -ne 0 ]; then
    cp $output_dir/mapping/${output_prefix}_avgMapping.prop.mat $BIN_DIR/final_warps_FS$fs_ver/
    cp $output_dir/mask/${template_sub_id}_cortex_estimate.nii.gz $BIN_DIR/liberal_cortex_masks_FS$fs_ver/
    echo "Warning: Average mapping & cortex mask have been copied to bin. Add -c 0 to command if you do not want to copy them."
  fi
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -s <template_sub_id> -n <num_of_sub> -r <RF_type> -l <ind_sub_list> -i <input_dir> -o <output_dir>

This script computes the average mapping from fsaverage to a volumetric atlas space, across a specified number of individual subjects, as step 3 in RF approaches. The averaging is carried out across the x/y/z index files previously projected from fsaverage to the volumetric atlas space through each individual subject.

After generating the average mapping, a loose cortex mask is also created using the count map generated in the process. This mask can be used in conjunction with the average mapping for final projection.

Note that both the average mapping and the cortex mask will be copied to the respective bin folders by default. Use -c option to specify otherwise.

REQUIRED ARGUMENTS:
	-s <template_sub_id> 	subject ID of the volumetric template used in recon-all

OPTIONAL ARGUMENTS:
	-n <num_of_sub>		number of subjects to use. This means taking the first <num_of_sub> subjects from <ind_sub_list>. For example, setting '-n 50' means the first 50 lines of <ind_sub_list> will be read to get subject IDs. Setting this to 0 will make the script use all subjects from <ind_sub_list>.
				[ default: 0 ]
	-r <RF_type>            type of RF approaches used in step 2, i.e. RF_M3Z or RF_ANTs
                                [ default: RF_ANTs ]
	-i <input_dir> 		absolute path to input directory. The inputs are the index files projected to the volumetric atlas space through individual subjects in step 2, i.e. the input directory should be the same as output directory in step 2.
				[ default: $(pwd)/results/index_\$template_sub_id ]
	-l <ind_sub_list> 	absolute path to a file containing individual subject IDs. Each line in the file should contain one subject ID.
				[ default: $DEFAULT_GSP_SUBLIST ]
	-o <output_dir> 	absolute path to output directory
				[ default: $(pwd)/results ]
	-c <copy_results>	set this to 0 to prevent the results generated from being copied to the bin folders. By default, both the average mapping and cortex mask generated will be copied to the respective bin folders.
				[ default: 1 ]
	-h			display help message

OUTPUTS:
	$0 will create 2 folders:
	1) mapping folder: 2 files will be generated, corresponding to the average mapping and count map in the volumetric atlas space. The count map shows at each votex how many subjects were projected to it.
	For example:
		allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_avgMapping.prop.mat
		allSub_fsaverage_to_FSL_MNI152_FS4.5.0_RF_ANTs_count.mat
	2) mask folder: 1 file will be generated, which is the liberal cortex mask generated using the count map
	For example:
		FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz

EXAMPLE:
	$0 -s 'FSL_MNI152_FS4.5.0'
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
RF_type=RF_ANTs
ind_sub_list=$DEFAULT_GSP_SUBLIST
output_dir=$(pwd)/results/
copy=1

#Assign arguments
while getopts "s:n:r:l:i:o:c:h" opt; do
  case $opt in
    s) template_sub_id=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    r) RF_type=${OPTARG} ;;
    l) ind_sub_list=${OPTARG} ;;
    i) input_dir=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    c) copy=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

if [ -z $input_dir ]; then
  input_dir=$(pwd)/results/index_$template_sub_id
fi

###########################################
#Check parameters
###########################################

if [ -z $template_sub_id ]; then
  echo "Template subject ID not defined."; 1>&2; exit 1
fi

###########################################
#Other set-ups
###########################################

#Make sure output directory is set up
if [ ! -d "$output_dir/mapping" ]; then
  echo "Output directory for mapping does not exist. Making directory now..."
  mkdir -p $output_dir/mapping
fi
if [ ! -d "$output_dir/mask" ]; then
  echo "Output directory for mask does not exist. Making directory now..."
  mkdir -p $output_dir/mask
fi

###########################################
#Implementation
###########################################

main



