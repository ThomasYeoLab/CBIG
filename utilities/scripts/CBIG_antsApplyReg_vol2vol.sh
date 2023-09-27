#!/usr/bin/env bash
# Written by Wu Jianxiao and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
# This function warps an input volumes to a target volume space with specified ANTs warp or inverse warp files

###########################################
# Main commands
###########################################
main(){

  # Set up warp files
  warp=${warp_dir}/${warp_prefix}1Warp.nii.gz
  inverse=${warp_dir}/${warp_prefix}1InverseWarp.nii.gz
  affine=${warp_dir}/${warp_prefix}0GenericAffine.mat

  # Set up fmri_reg command
  if [ ! -z $fmri_reg ]; then
    fmri_warp="-t $fmri_reg"
  fi

  # Set up time_series command
  if [ $time_series -eq 1 ]; then
    vol4d="-e 3"
  fi

  # Warp input volume to target's T1 space
  output=$output_dir/${output_prefix}.nii.gz
  if [ ! -e $output ]; then
    case $warp_setting in
      forward)
        cmd="${ANTs_dir}/antsApplyTransforms -d 3 $vol4d -i $input -r $target -n $interp -t $warp -t $affine $fmri_warp"
        cmd="$cmd -o $output"
        echo $cmd
        eval $cmd
        ;;
      inverse)
        cmd="${ANTs_dir}/antsApplyTransforms -d 3 $vol4d -i $input -r $target -n $interp $fmri_warp -t [$affine, 1]"
        cmd="$cmd -t $inverse -o $output"
        echo $cmd
        eval $cmd
        ;;
      *)
        echo "Invalid warp setting. Use either 'forward' or 'inverse'"
     esac
  else
    echo "The input volume have already been projected into ${target}'s space by ANTs"
  fi
}

##################################################################
# Function usage
##################################################################

# Usage
usage() { echo "
Usage: CBIG_antsApplyReg_vol2vol.sh -i input -r target -w warp_prefix -p output_prefix

This script applies existing warps prepared using ANTs registration to map an input volume to a target volume.

REQUIRED ARGUMENTS:

    -i <input>           absolute path to input volume to map from
    -r <target>          absolute path to target volume to be mapped to
    -w <warp_prefix>     prefix of the warp files. Specifically, the affine, forward and inverse warps should have the 
                         following names:
                             warp_prefix0GenericAffine.mat
                             warp_prefix1Warp.nii.gz
                             warp_prefix1InverseWarp.nii.gz
    -p <output_prefix>   prefix for output volume

OPTIONAL ARGUMENTS:
    -f <fmri_reg>        absolute path to native-to-T1 registration file for fMRI input. The file must have been 
                         converted to format used by ANTs (e.g. .txt or .mat) if it was generated using other tools.
                         [ default: unset ]
    -e <time_series>     set this to 1 if the input volume is 4D (time-series)
                         [ default: 0 ]
    -d <warp_dir>        absolute path to the warps
                         [ default: $(pwd)/results ]
    -o <output_dir>      absolute path to output directory
                         [ default: $(pwd)/results ]
    -s <warp setting>    warping direction ('forward' for forward warp, 'inverse' for inverse warp). For example, if 
                         the target space was registered to the input space during ANTs registration, then inverse 
                         warp should be used to now map the input to the target.
                         [ default: forward ]
    -t <interp>          interpolation (Linear, NearestNeighbor, etc.)
                         [ default: Linear ]
    -a <ANTs_dir>        directory where ANTs is installed
                         [ default: $CBIG_ANTS_DIR ]
    -h                   display help message

OUTPUTS:
    $0 will create 1 output file in the output directory, corresponding to the input warped into target's space.
    For example:
        output_prefix.nii.gz

EXAMPLE:
    $0 -i /path/to/my/data.nii.gz -r /path/to/my/template.nii.gz -w subject_moving_template_fixed -p data_to_template 
    -f data_to_T1.txt -s forward

" 1>&2; exit 1; }

# Display help message if no argument is supplied
if [ $# -eq 0 ]; then
  usage; 1>&2; exit 1
fi

##################################################################
# Assign input variables
##################################################################

# Default parameter
warp_dir=$(pwd)/results
output_dir=$(pwd)/results
warp_setting=forward
interp=Linear
ANTs_dir=$CBIG_ANTS_DIR
time_series=0

# Assign parameter
while getopts "i:r:w:p:f:e:d:o:s:t:a:h" opt; do
  case $opt in
    i) input=${OPTARG} ;;
    r) target=${OPTARG} ;;
    w) warp_prefix=${OPTARG} ;;
    p) output_prefix=${OPTARG} ;;
    f) fmri_reg=${OPTARG} ;;
    e) time_series=${OPTARG} ;;
    d) warp_dir=${OPTARG} ;;
    o) output_dir=${OPTARG} ;;
    s) warp_setting=${OPTARG} ;;
    t) interp=${OPTARG} ;;
    a) ANTs_dir=${OPTARG} ;;
    h) usage; exit ;;
    *) usage; 1>&2; exit 1 ;;
  esac
done

##################################################################
# Check parameter
##################################################################

if [ -z $target ]; then
  echo "Reference volume not defined."; 1>&2; exit 1
fi

if [ -z $input ]; then
  echo "Input volume not defined."; 1>&2; exit 1
fi

if [ -z $warp_prefix ]; then
  echo "Warp prefix not defined."; 1>&2; exit 1
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


