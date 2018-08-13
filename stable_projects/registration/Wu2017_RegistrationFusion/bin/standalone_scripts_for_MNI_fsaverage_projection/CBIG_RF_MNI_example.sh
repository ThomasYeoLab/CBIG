#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This warpper script runs CBIG_projectVol2fsaverage.sh for an example input

if [ "$(uname)" == "Linux" ]; then
  SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
elif [ "$(uname)" == "Darwin" ]; then
  SCRIPT_DIR=$(dirname "$0")
  echo "Note: path to scripts may not be retrieved properly sometimes (e.g. when symlinks are involved) on Mac. Either make sure you are not calling a symlink version of the script, or call the script from the stand-alone folder itself."
fi

###########################################
# Main commands
###########################################
main(){
  # Set-ups
  input_name=MNI_probMap_ants.central_sulc
  input=$SCRIPT_DIR/${input_name}.nii.gz

  # Call vol2fsaverage script
  cmd="$SCRIPT_DIR/CBIG_RF_projectMNI2fsaverage.sh -s $input -f $FREESURFER_HOME -m $matlab_bin -o $(pwd)"
  echo $cmd
  eval $cmd

  # Visualise lh results
  output=lh.${input_name}.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz
  demo=$SCRIPT_DIR/lh.${input_name}.demo.nii.gz
  cmd="freeview -f $FREESURFER_HOME/subjects/fsaverage/surf/lh.inflated:overlay=$output:overlay_threshold=0.01,1:overlay=$demo:overlay_threshold=0.01,1 -viewport '3d'"
  echo $cmd
  eval $cmd
}

###########################################
# Function usage
###########################################

# usage
usage() { echo "
Usage: $0 -f <freesurfer_dir> -m <matlab_bin>

This script projects an example input in MNI152 space to fsaverage space and visualise the results in freeview. 

OPTIONAL ARGUMENTS:
	-f <freesurfer_dir>	absolute path to FreeSurfer directory (this will overwrite the default \$FREESUFER_HOME). This setting may affect the version of FreeSurfer .m scripts called, but not the version of the warps/mappings used.
				[ default: $FREESURFER_HOME ]
	-m <matlab_bin>		absolute path to the Matlab bin directory. Default is set to \$CBIG_MATLAB_DIR/bin or empty. This setting may affect the version of Matlab builtin .m scripts called.
				[ default: $CBIG_MATLAB_DIR/bin ]
	-h			display help message

OUTPUTS:
	$0 will create 2 files, corresponding to the projected data onto fsaverage left and right hemispheres respectively, i.e.:
		lh.MNI_probMap_ants.central_sulc.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz
		rh.MNI_probMap_ants.central_sulc.allSub_RF_ANTs_MNI152_orig_to_fsaverage.nii.gz

" 1>&2; exit 1; }

###########################################
# Parse arguments
###########################################

matlab_bin=$CBIG_MATLAB_DIR/bin
while getopts "f:m:h" opt; do
  case $opt in
    f) export FREESURFER_HOME=${OPTARG} ;;
    m) matlab_bin=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

###########################################
# Check parameters
###########################################

if [ ! -d $FREESURFER_HOME ]; then
  echo "FreeSurfer directory cannot be found. You can use -f to manually set the path to FreeSurfer if necessary."; 1>&2; exit 1
fi

if [ ! -e $matlab_bin/matlab ]; then
  echo "Matlab cannot be found. You can use -m to manually set the path to Matlab bin if necessary."; 1>&2; exit 1
fi

###########################################
# Implementation
###########################################

main
