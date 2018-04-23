#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#This script uses the mappings generate using Registration Fusion approach to project batch data in a volumetric space to fsaverage surface

###########################################
#Set-up
###########################################

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
fs_ver=`cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' | head -c 3`
WARP_DIR=$(dirname "$SCRIPT_DIR")/final_warps_FS$fs_ver
MAKE_VOLINDEX_SCRIPT=$(dirname "$(dirname "$SCRIPT_DIR")")/registration_fusion/scripts_vol2surf/CBIG_RF_step1_make_xyzIndex_volTemplate.sh

###########################################
#Main commands
###########################################

main(){
  #Set up file prefix
  if [ $num_sub -gt 0 ]; then
    prefix=${num_sub}Sub_${RF_type}_${template_type}_to_fsaverage
    else
      prefix=allSub_${RF_type}_${template_type}_to_fsaverage
  fi

  #Set up file names
  lh_map=$map_dir/lh.avgMapping_$prefix.mat
  rh_map=$map_dir/rh.avgMapping_$prefix.mat
  input_paths=`cat $input_list`

  #Loop through each input
  for input in $input_paths
  do
    #Set up output file name
    input_name=`basename $input .gz`
    input_name=`basename $input_name .nii`
    lh_output=$output_dir/lh.$input_name.$prefix.nii.gz
    rh_output=$output_dir/rh.$input_name.$prefix.nii.gz

    #Call function to run projecting of input
    matlab -nodesktop -nojvm -nosplash -r "mri = MRIread('$input'); \
                                           addpath('$SCRIPT_DIR'); \
                                           [lh_proj_data, rh_proj_data] = CBIG_RF_projectVol2fsaverage_single(mri, '$interp', '$lh_map', '$rh_map', '$average'); \
                                           mri.vol = permute(lh_proj_data, [4 2 3 1]); \
                                           MRIwrite(mri, '$lh_output'); \
                                           mri.vol = permute(rh_proj_data, [4 2 3 1]); \
                                           MRIwrite(mri, '$rh_output'); \
                                           exit"
  done
}

###########################################
#Function usage
###########################################

#usage
usage() { echo "
Usage: $0 -l <input_list> -p <template_type> -n <num_of_sub> -r <RF_type> -d <map_dir> -i <interp> -f <average_mesh> -o <output_dir>

This script projects a batch of new data in a volumetric atlas space to fsaverage space, using the mapping generated using Registration Fusion approach. Note that the default version of mapping used is determined by the current FreeSurfer version used.

REQUIRED ARGUMENTS:
	-l <input_list>         absolute path to file containing input data file names. Each line in the file should contain full path to one file. Note that inputs should be in the format of .nii or .nii.gz

OPTIONAL ARGUMENTS:
	-p <template_type>      type of volumetric template used in index files creation. See $MAKE_VOLINDEX_SCRIPT for more details.
				[ default: MNI152_orig ]
	-n <num_of_sub>		number of subjects used in creating the average mapping. For example, setting '-n 50' means the mapping used was averaged across 50 subjects. Setting this to 0 will make the script use the mapping averaged across all subjects.
				[ default: 0 ]
	-r <RF_type>            type of Registration Fusion approaches used to generate the mappings (RF_M3Z or RF_ANTs). RF-M3Z is recommended if data was registered from subject's space to the volumetric atlas space using FreeSurfer. RF-ANTs is recommended if such registrations were carried out using other tools, especially ANTs.
				[ default: RF_ANTs ]
	-d <map_dir>          	absolute path to mapping directory. The mappings are the average mappings generated in Registration Fusion approaches. The default average mapping used all GSP subjects with FreeSurfer $fs_ver. 
				[ default: $WARP_DIR ]
	-i <interp> 		interpolation (linear or nearest)
				[ default: linear ]
	-f <average_mesh>	fsaverage mesh version (fsaverage, fsaverage5 or fsaverage6)
				[ default: fsaverage ]
	-o <output_dir> 	absolute path to output directory
				[ default: $(pwd)/results/projected_vol2fsaverage ]
	-h			display help message

OUTPUTS:
	$0 will create 2 files for each input, corresponding to the projected data onto fsaverage left and right hemispheres respectively. 
	For example: 
		lh.input_name.allSub_RF_ANTs_MNI152_norm_to_fsaverage.nii.gz
		rh.input_name.allSub_RF_ANTs_MNI152_norm_to_fsaverage.nii.gz

EXAMPLE:
	$0 -l my_input_list.csv
	$0 -l my_input_list.csv -p Colin27_norm -r RF_M3Z -n 50

" 1>&2; exit 1; }

#Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

###########################################
#Parse arguments
###########################################

#Default parameters
template_type=MNI152_orig
num_sub=0
RF_type=RF_ANTs
map_dir=$WARP_DIR
interp=linear
average=fsaverage
output_dir=$(pwd)/results/projected_vol2fsaverage/

#Assign arguments
while getopts "l:p:n:r:d:i:f:o:h" opt; do
  case $opt in
    l) input_list=${OPTARG} ;;
    p) template_type=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    r) RF_type=${OPTARG} ;;
    d) map_dir=${OPTARG} ;;
    i) interp=${OPTARG} ;;
    f) average=${OPTARG} ;;
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
