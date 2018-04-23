#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This function uses the mapping generated using Registration Fusion approach to project batch surface data in fsaverage to a volumetric atlas space

###########################################
#Define paths
###########################################

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
fs_ver=`cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' | head -c 3`
WARP_DIR=$(dirname "$SCRIPT_DIR")/final_warps_FS$fs_ver
DEFAULT_MASK=$(dirname "$SCRIPT_DIR")/liberal_cortex_masks_FS$fs_ver/FSL_MNI152_FS4.5.0_cortex_estimate.nii.gz

###########################################
#Main commands
###########################################

main(){
  #Set up file prefix
  if [ $num_sub -gt 0 ]; then
    prefix=${num_sub}Sub_fsaverage_to_${template_sub_id}_$RF_type
    else
      prefix=allSub_fsaverage_to_${template_sub_id}_$RF_type
  fi

  #Set up file names
  input_paths=`cat $input_list`
  map=${map_dir}/${prefix}_avgMapping.prop.mat

  #Loop through each input
  for input in $input_paths
  do
    #Set up output file name
    input_name=`basename $input .mat`
    output=$output_dir/$input_name.$prefix.nii.gz
    output_seg=$output_dir/seg.$input_name.$prefix.nii.gz

    #Call function to run projecting of input
    matlab -nojvm -nosplash -nodesktop -r "load('$input'); \
                                           addpath('$SCRIPT_DIR'); \
                                           [projected, projected_seg] = CBIG_RF_projectfsaverage2Vol_single(lh_label, rh_label, '$interp', '$map', '$mask'); \
                                           MRIwrite(projected, '$output'); \
                                           MRIwrite(projected_seg, '$output_seg'); \
                                           exit"
  done
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -n <num_of_sub> -l <input_list> -d <map_dir> -p <avg_map_prefix> -i <interp> -m <cortical_mask> -o <output_dir>

This script projects a batch of new data in fsaverage to a volumetric atlas space, using the mapping generated using Registration Fusion approach. Note that the default version of mapping and mask used is determined by the current FreeSurfer version used.

REQUIRED ARGUMENTS:
	-l <input_list>         absolute path to file containing input data file names. Each line in the file should contain full path to one file. Note that the inputs should be in .mat format, containing two variables 'lh_label' and 'rh_label', corresponding to the input data in left and right hemispheres respectively.

OPTIONAL ARGUMENTS:
	-s <template_sub_id> 	subject ID of the volumetric template used in recon-all
				[ default: FSL_MNI152_FS4.5.0 ]
	-n <num_of_sub>		number of subjects used in creating the average mapping. For example, setting '-n 50' means the mapping used was averaged across 50 subjects. Setting this to 0 will make the script use the mapping averaged across all subjects.
				[ default: 0 ]
	-r <RF_type>            type of Registration Fusion approaches used to generate the mapping (RF_M3Z or RF_ANTs). RF-M3Z is recommended if data was registered from subject's space to the volumetric atlas space using FreeSurfer. RF-ANTs is recommended if such registrations were carried out using other tools, especially ANTs.
				[ default: RF_ANTs ]
	-d <map_dir> 		absolute path to mapping directory. The mappings are the average mappings generated in Registration Fusion approaches. The default average mapping used all GSP subjects with FreeSurfer $fs_ver. 
				[ default: $WARP_DIR ]
	-i <interp> 		interpolation (linear or nearest)
				[ default: nearest ]
	-m <cortical_mask> 	absolute path to liberal cortical mask. This should be a binary mask that loosely defines the cortex of the volumetric atlas space. The shape of the mask will determine the shape of the output. The default mask is generated in Registration Fusion approach with FreeSurfer $fs_ver.
				[ default: $DEFAULT_MASK ]
	-o <output_dir> 	absolute path to output directory
				[ default: $(pwd)/results/projected_fsaverage2vol ]
	-h			display help message

OUTPUTS:
	$0 will create 2 files for each input, corresponding to the projected data in the volumetric atlas sapce and the projected data in segmentation form (with left hemisphere values starting from 0 and right hemisphere values starting from 1000). 
	For example: 
		input_name.allSub_fsaverage_to_FSL_MNI152_FS4.5_RF_ANTs.nii.gz
		seg.input_name.allSub_fsaverage_to_FSL_MNI152_FS4.5_RF_ANTs.nii.gz

EXAMPLE:
	$0 -l my_input_list.csv
	$0 -l my_input_list.csv -s SPM_Colin27_FS4.5.0 -n 50 

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

###########################################
#Parse arguments
###########################################

#Default parameters
template_sub_id=FSL_MNI152_FS4.5.0
num_sub=0
RF_type=RF_ANTs
map_dir=$WARP_DIR
interp=nearest
mask=$DEFAULT_MASK
output_dir=$(pwd)/results/projected_fsaverage2vol

#Assign arguments
while getopts "l:s:n:r:d:i:m:o:h" opt; do
  case $opt in
    s) template_sub_id=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    l) input_list=${OPTARG} ;;
    r) RF_type=${OPTARG} ;;
    d) map_dir=${OPTARG} ;;
    i) interp=${OPTARG} ;;
    m) mask=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

###########################################
#Check parameters
###########################################

if [ -z $input_list ]; then
  echo "Input list not defined."; 1>&2; exit 1
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



