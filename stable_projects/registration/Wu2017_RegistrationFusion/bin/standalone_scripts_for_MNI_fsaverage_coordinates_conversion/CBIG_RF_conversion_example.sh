#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This warpper script runs CBIG_RF_MNICoord2fsaverageVertex.m and CBIG_RF_fsaverageVertex2MNICoord.m for an example

if [ "$(uname)" == "Linux" ]; then
  SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
elif [ "$(uname)" == "Darwin" ]; then
  SCRIPT_DIR=$(dirname "$0")
  echo "Note: path to scripts may not be retrieved properly sometimes (e.g. when symlinks are involved) on Mac. \\
  Either make sure you are not calling a symlink version of the script, or call the script from the stand-alone \\
  folder itself."
fi

###########################################
# Main commands
###########################################
main(){

  # Call the conversion functions in Matlab
  $matlab_bin/matlab -nodesktop -nojvm -nosplash -r "mni_coords = [60; 0; 10]; \
    fprintf('Starting set of MNI152 coordinates: [%d, %d, %d]\n', mni_coords); \
    vertices = CBIG_RF_MNICoord2fsaverageVertex(mni_coords); \
    fprintf('The corresponding vertex number is: %d\n', vertices); \
    mni_coords_back = CBIG_RF_fsaverageVertex2MNICoord('rh', vertices); \
    fprintf('MNI152 coordinates corresponding to vertex %d are: [%.2f, %.2f, %.2f]\n', vertices, mni_coords_back); \
    dist = sqrt(sum((mni_coords - mni_coords_back).^2)); \
    fprintf('The distance between original and converted-back MNI152 coordinates is %.2f mm \n', dist); \
    exit"

}

###########################################
# Function usage
###########################################

# usage
usage() { echo "
Usage: $0 -f <freesurfer_dir> -m <matlab_bin>

This script converts an example set of RAS coordinates ([60, 0, 10]) in MNI152 space to the corresponding 
fsaverage vertex, and back to coordinates in MNI152 space. 

OPTIONAL ARGUMENTS:
	-f <freesurfer_dir>	absolute path to FreeSurfer directory (this will overwrite the default \$FREESUFER_HOME).
				[ default: $FREESURFER_HOME ]
	-m <matlab_bin>		absolute path to the Matlab bin directory. Default is set to \$CBIG_MATLAB_DIR/bin or empty. 
                    This setting may affect the version of Matlab builtin .m scripts called.
				[ default: $CBIG_MATLAB_DIR/bin ]
	-h			display help message

OUTPUTS:
	$0 does not create any output file. The results will be printed on screen.

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
  echo "FreeSurfer directory cannot be found. You can use -f to manually set the path to FreeSurfer if necessary." 
  1>&2; exit 1
fi

if [ ! -e $matlab_bin/matlab ]; then
  echo "Matlab cannot be found. You can use -m to manually set the path to Matlab bin if necessary."; 1>&2; exit 1
fi

###########################################
# Implementation
###########################################

main
