#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This script uses the mappings generate using Registration Fusion approach to project batch data in a volumetric space to fsaverage surface

###########################################
# Define paths
###########################################

if [ "$(uname)" == "Linux" ]; then
  SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
elif [ "$(uname)" == "Darwin" ]; then
  SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
  echo "Note: path to scripts may not be retrieved properly sometimes (e.g. when symlinks are involved) on Mac. Either make sure you are not calling a symlink version of the script, or call the script from the stand-alone folder itself."
fi
WARP_DIR=$(dirname "$SCRIPT_DIR")/final_warps_FS5.3

###########################################
# Main commands
###########################################

main(){
  # Set up file prefix
  if [ $num_sub -gt 0 ]; then
    prefix=${num_sub}Sub_${RF_type}_${template_type}_to_fsaverage
    else
      prefix=allSub_${RF_type}_${template_type}_to_fsaverage
  fi

  # Set up file names
  lh_map=$map_dir/lh.avgMapping_$prefix.mat
  rh_map=$map_dir/rh.avgMapping_$prefix.mat

  # Set up output file name
  input_name=`basename $input .gz`
  input_name=`basename $input_name .nii`
  lh_output=$output_dir/lh.$input_name.$prefix.nii.gz
  rh_output=$output_dir/rh.$input_name.$prefix.nii.gz

  # Call function to run projecting of input
  $matlab_bin/matlab -nodesktop -nojvm -nosplash -r "addpath('$SCRIPT_DIR'); \
                                         input = MRIread('$input'); \
                                         [lh_proj_data, rh_proj_data] = CBIG_RF_projectMNI2fsaverage('$input', '$interp', '$lh_map', '$rh_map'); \
                                         input.vol = permute(lh_proj_data, [4 2 3 1]); \
                                         MRIwrite(input, '$lh_output'); \
                                         input.vol = permute(rh_proj_data, [4 2 3 1]); \
                                         MRIwrite(input, '$rh_output'); \
                                         exit"
}

###########################################
# Function usage
###########################################

# usage
usage() { echo "
Usage: $0 -s <input> -o <output_dir> -p <template_type> -n <num_of_sub> -r <RF_type> -d <map_dir> -i <interp> -f <freesurfer_dir> -m <matlab_bin>

This script projects data in a volumetric atlas space to fsaverage space, using the mapping generated using Registration Fusion approach. 

REQUIRED ARGUMENTS:
	-s <input>              absolute path to input volume. Input should be in nifti format.
	-o <output_dir>         absolute path to output directory.

OPTIONAL ARGUMENTS:
	-p <template_type>      type of volumetric template used in index files creation. Use 'MNI152_norm' or 'Colin27_norm' for RF-M3Z MNI152-to-fsaverage or Colin27-to-fsaverage mappings, 'MNI152_orig' or 'Colin27_orig' for the respective RF-ANTs mappings.
	                        [ default: MNI152_orig ]
	-n <num_of_sub>         number of subjects used in creating the average mapping. For example, setting '-n 50' means the mapping used was averaged across 50 subjects. Setting this to 0 will make the script use the mapping averaged across all subjects.
	                        [ default: 0 ]
	-r <RF_type>            type of Registration Fusion approaches used to generate the mappings (RF_M3Z or RF_ANTs). RF-M3Z is recommended if data was registered from subject's space to the volumetric atlas space using FreeSurfer. RF-ANTs is recommended if such registrations were carried out using other tools, especially ANTs.
	                        [ default: RF_ANTs ]
	-d <map_dir>            absolute path to mapping directory. The mappings are the average mappings generated in Registration Fusion approaches. The average mapping using all GSP subjects are used by default. 
	                        [ default: $WARP_DIR ]
	-i <interp>             interpolation (linear or nearest)
	                        [ default: linear ]
	-f <freesurfer_dir>     absolute path to FreeSurfer directory (this will overwrite the default \$FREESUFER_HOME). This setting may affect the version of FreeSurfer .m scripts called, but not the version of the warps/mappings used.
	                        [ default: $FREESURFER_HOME ]
	-m <matlab_bin>         absolute path to the Matlab bin directory. Default is set to \$CBIG_MATLAB_DIR/bin or empty. This setting may affect the version of Matlab builtin .m scripts called.
	                        [ default: $CBIG_MATLAB_DIR/bin ]
	-h                      display help message

OUTPUTS:
	$0 will create 2 files, corresponding to the projected data onto fsaverage left and right hemispheres respectively. 
	For example: 
		lh.input_name.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz
		rh.input_name.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz

EXAMPLE:
	$0 -s $SCRIPT_DIR/MNI_probMap_ants.central_sulc.nii.gz -o \$(pwd)
	freeview -f \$FREESURFER_HOME/subjects/fsaverage/surf/lh.inflated:overlay=lh.MNI_probMap_ants.central_sulc.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz:overlay=${SCRIPT_DIR}/lh.MNI_probMap_ants.central_sulc.demo.nii.gz:overlay_threshold=0.001,1 -viewport '3d'

	This projects the MNI probabilistic map of central sulcus to fsaverage surface and visualise the left hemisphere results in freeview.

" 1>&2; exit 1; }

# Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

###########################################
# Parse arguments
###########################################

# Default parameters
template_type=MNI152_orig
num_sub=0
RF_type=RF_ANTs
map_dir=$WARP_DIR
interp=linear
matlab_bin=$CBIG_MATLAB_DIR/bin

# Assign arguments
while getopts "s:p:n:r:d:i:o:f:m:h" opt; do
  case $opt in
    s) input=${OPTARG} ;;
    p) template_type=${OPTARG} ;;
    n) num_sub=${OPTARG} ;;
    r) RF_type=${OPTARG} ;;
    d) map_dir=${OPTARG} ;;
    i) interp=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    f) export FREESURFER_HOME=${OPTARG} ;;
    m) matlab_bin=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

###########################################
# Check parameters
###########################################

if [ -z $input ]; then
  echo "Input not defined."; 1>&2; exit 1
fi

if [ -z $output_dir ]; then
  echo "Output directory not defined."; 1>&2; exit 1
fi

if [ ! -d $FREESURFER_HOME ]; then
  echo "FreeSurfer directory cannot be found. You can use -f to manually set the path to FreeSurfer if necessary."; 1>&2; exit 1
else
  fs_ver=`cat $FREESURFER_HOME/build-stamp.txt | sed 's@.*-v@@' | sed 's@-.*@@' | head -c 3`
  echo "Using FreeSurfer version $fs_ver"
fi

if [ ! -e $matlab_bin/matlab ]; then
  echo "Matlab cannot be found. You can use -m to manually set the path to Matlab bin if necessary."; 1>&2; exit 1
fi

###########################################
# Other set-ups
###########################################

# Make sure output directory is set up
if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
# Implementation
###########################################

main
