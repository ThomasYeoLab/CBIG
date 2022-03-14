#!/bin/sh
# This script generate a "super inflated" surface from very_inflated surface. 
# This surface is more inflated than the very_inflated version, which cannot 
# show the insula very well. You can load the output files in wb_view.
# See https://www.humanconnectome.org/software/workbench-command/-surface-inflation
# for more details.
#
# Input:
#       output_dir: Path of output folder
#
# Output:
#       $output_dir/fsaverage.L.super_inflated.32k_fs_LR.surf.gii
#       $output_dir/fsaverage.R.super_inflated.32k_fs_LR.surf.gii
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

output_dir=$1
input_dir="$CBIG_CODE_DIR/data/templates/surface/fs_LR_32k"
lh_input_surf="${input_dir}/fsaverage.L.very_inflated.32k_fs_LR.surf.gii"
rh_input_surf="${input_dir}/fsaverage.R.very_inflated.32k_fs_LR.surf.gii"
lh_ana_surf="${input_dir}/fsaverage.L.pial_orig.32k_fs_LR.surf.gii"
rh_ana_surf="${input_dir}/fsaverage.R.pial_orig.32k_fs_LR.surf.gii"

lh_output_surf="${output_dir}/fsaverage.L.super_inflated.32k_fs_LR.surf.gii"
rh_output_surf="${output_dir}/fsaverage.R.super_inflated.32k_fs_LR.surf.gii"
wb_command -surface-inflation ${lh_ana_surf} ${lh_input_surf} 10 0.9 20 1.02 ${lh_output_surf}
wb_command -surface-inflation ${rh_ana_surf} ${rh_input_surf} 10 0.9 20 1.02 ${rh_output_surf}
