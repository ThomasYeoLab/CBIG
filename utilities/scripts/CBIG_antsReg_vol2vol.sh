#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This function registers a moving volume to a fixed volume using ANTs

###########################################
# Main commands
###########################################
main(){
  # Refer to the official example (https://github.com/stnava/ANTs/blob/master/Scripts/newAntsExample.sh) for choice of 
  # parameters. 
  # In this script, compared to the example production setting, the translation and rigid stages are removed, and max number 
  # of iterations at affine stage reduced.

  output=$output_dir/${output_prefix} 
  cmd="${ANTs_dir}/antsRegistration -d 3 -r [ $fixed, $moving, 1 ] \
      -m mattes[ $fixed, $moving, 1, 32, regular, 0.3] \
        -t affine[ 0.1 ] \
        -c [ $iter_affine, 1.e-8, 20 ] \
        -s 4x2x1vox \
        -f 3x2x1 \
      -m cc[ $fixed, $moving, 1, 4] \
        -t SyN[ .20, 3, 0] \
        -c [ $iter_SyN, 0, 5 ] \
        -s 1x0.5x0vox \
        -f 4x2x1 -u 1 -z 1\
      -o ${output}"
   echo $cmd
   eval $cmd
}

##################################################################
# Function usage
##################################################################

# Usage
usage() { echo "
Usage: CBIG_antsReg_vol2vol.sh -r reference -i input -p output_prefix

This script runs ANTs registration from an input volume to a reference volume, generating a forward warp, an inverse warp and
an affine transform files.

REQUIRED ARGUMENTS:
    -r <reference>        absolute path to volume to be mapped to
    -i <input>            absolute path to volume to map from
    -p <output_prefix>    prefix for output warps

OPTIONAL ARGUMENTS:
    -d <output_dir>       absolute path to output directory
                          [ default: $(pwd)/results ]
    -a <ANTs_dir>         directory where ANTs is installed 
                          [ default: $CBIG_ANTS_DIR ]
    -j <iterations_of_affine_transform>
                          [ default: 100x100x200 ]
    -k <iterations_of_SyN_transform>
                          [ default: 100x100x50 ]
    -h                    display help message

OUTPUTS:
    $0 will create 3 warp files in the output directory, corresponding to the affine, forward and inverse warps.
    These outputs can be supplied to $CBIG_CODE_DIR/utilities/scripts/CBIG_antsApplyReg_vol2vol.sh to warp data in the 
input space to the reference space.
    For example:
        output_prefix0GenericAffine.mat
        output_prefix1Warp.nii.gz
        output_prefix1InverseWarp.nii.gz

EXAMPLE:
    $0 -r /path/to/my/template.nii.gz -i /path/to/my/subject.nii.gz -p subject_moving_template_fixed

" 1>&2; exit 1; }

# Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
# Assign input variables
##################################################################

# Default parameter
output_dir=$(pwd)/results
ANTs_dir=$CBIG_ANTS_DIR
iter_affine=100x100x200
iter_SyN=100x100x50

# Assign parameter
while getopts "r:i:d:p:a:h:j:k:" opt; do
  case $opt in
    r) fixed=${OPTARG} ;;
    i) moving=${OPTARG} ;;
    d) output_dir=${OPTARG} ;;
    p) output_prefix=${OPTARG} ;;
    a) ANTs_dir=${OPTARG} ;;
    j) iter_affine=${OPTARG} ;;
    k) iter_SyN=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

echo "iter_affine = $iter_affine"
echo "iter_SyN = $iter_SyN"

##################################################################
# Check parameter
##################################################################

if [ -z $fixed ]; then
  echo "Reference volume not defined."; 1>&2; exit 1
fi

if [ -z $moving ]; then
  echo "Input volume not defined."; 1>&2; exit 1
fi

if [ -z $output_prefix ]; then
  echo "Output prefix not defined."; 1>&2; exit 1
fi

##################################################################
# Disable multi-threading
##################################################################

ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS

##################################################################
# Set up output directory
##################################################################

if [ ! -d "$output_dir" ]; then
  echo "Output directory does not exist. Making directory now..."
  mkdir -p $output_dir
fi

###########################################
# Implementation
###########################################

main



