#!/bin/sh
# CBIG_IndCBM_create_template.sh <surf_mesh> <binary_mask> <output_dir> -m(optional) <midline> \
# -l(optional) <lh_vertex> -r(optional) <rh_vertex> 
# This function creates a cifti template combining the cerebral cortical surface and the cerebellum in the volume. 
# Input:
#   surf_mesh:      Mesh neame for your surface. Can be fsaverage5, fsaverage6, fsaverage, fs_LR_32k, fs_LR_164k or
#                   other customized mesh. Need to pass in optional input '-l' and '-r' when using customized mesh
#                   name. 
#   binary_mask:    Mask of the cerebellum in the volume (nii or nii.gz). Cerebellum should be labeled as 1 and rest 
#                   should be 0.
#   output_dir:     Path of output folder. 
# Optional Input:
#   -m <midline>:   This input is optional. Default setting assume the left hemisphere and right hemisphere of the 
#                   cerebellum are divided at the middle of the volume. You can specify x=? as the dividing line. 
#                   Shift in this only causes different voxel ordering in the template file and does not affect the 
#                   visualization of the final parcellation output.
#   -l <lh_vertex>: Number of vertices of the left hemisphere. Only need to pass in when mesh name is not fsaverage5, 
#                   fsaverage6, fsaverage, fs_LR_32k or fs_LR_164k.
#   -r <lh_vertex>: Number of vertices of the right hemisphere. Only need to pass in when mesh name is not fsaverage5, 
#                   fsaverage6, fsaverage, fs_LR_32k or fs_LR_164k.
# Output:
#                   Generated template file will be saved under the output folder as 
#                   <binary_mask_file_name>_<surf_mesh>_cerebellum_template.dscalar.nii
# Example:  CBIG_IndCBM_create_template.sh fsaverage5 <binary_mask> <output_dir> -m 90
#           CBIG_IndCBM_create_template.sh individual_surf <binary_mask> <output_dir> -l 123456 -r 123457
# Written by XUE Aihuiping and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

surf_mesh=$1
bin_mask_file=$2
output_dir=$3

# Handle optional inputs 
OPTIND=$(($OPTIND+3))
while getopts "m:l:r:" opt; do
    case $opt in
        m) midline=${OPTARG};;
        l) N_lh=${OPTARG};;
        r) N_rh=${OPTARG};;
        ?) echo "Unknown parameter";;
    esac
done

if [ -n "$N_lh" ];then
    echo "Left hemispere has ${N_lh} vertices."
fi
if [ -n "$N_rh" ];then
    echo "Right hemispere has ${N_rh} vertices."
fi

# Compute midline by picking the middle slice
if [ ! -n "$midline" ] ;then
    midline=`mri_info ${bin_mask_file} --ncols`
    midline=$((midline/2-1))
    echo "Midline automatically set at x=$midline."
fi

if [ ! -d ${output_dir} ];then
    mkdir -p ${output_dir}
fi

name=`basename ${bin_mask_file} .nii.gz`
script_dir=`dirname $(readlink -f $0)`
LUT="${script_dir}/lib/CBIG_IndCBM_cerebellum_mask_LUT.txt"
template="${output_dir}/${name}_${surf_mesh}_cerebellum_template.dscalar.nii"

echo "Creating template for ${name}."

lh_file="${output_dir}/${name}_lh_mask.nii.gz"
if [ -f ${lh_file} ];then
    rm -f ${lh_file}
fi
fslmaths ${bin_mask_file} -roi ${midline} -1 0 -1 0 -1 0 -1 ${lh_file}

rh_file="${output_dir}/${name}_rh_mask.nii.gz"
if [ -f ${rh_file} ];then
    rm -f ${rh_file}
fi

# Left hemisphere of cerebellum labeled as 8 and right hemisphere labeled as 47
fslmaths ${bin_mask_file} -roi 0 ${midline} 0 -1 0 -1 0 -1 ${rh_file}
wb_command -volume-math '8*x' ${lh_file} -var x ${lh_file}
wb_command -volume-math '47*x' ${rh_file} -var x ${rh_file}

mask_file="${output_dir}/${name}_mask.nii.gz"
if [ -f ${mask_file} ];then
    rm -f ${mask_file}
fi
wb_command -volume-math 'x + y' ${mask_file} -var x ${lh_file} -var y ${rh_file}

label_mask_file="${output_dir}/${name}_label_mask.nii.gz"
if [ -f ${label_mask_file} ];then
    rm -f ${label_mask_file}
fi
wb_command -volume-label-import ${mask_file} ${LUT} ${label_mask_file} -discard-others

# Generate surface gii files
if [ -n "$N_lh" ] && [ -n "$N_rh" ];then
    matlab -nosplash -nodisplay -nodesktop -r " \
CBIG_CODE_DIR=getenv('CBIG_CODE_DIR'); \
addpath(fullfile(CBIG_CODE_DIR, 'external_packages', 'matlab', 'default_packages', 'cifti-matlab')); \
addpath(genpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Xue2021_IndCerebellum', 'lib')));\
CBIG_IndCBM_create_surf_gifti('${surf_mesh}', '${script_dir}', 'lh', '${N_lh}', 'rh', '${N_rh}'); \
exit;"
else
    matlab -nosplash -nodisplay -nodesktop -r " \
CBIG_CODE_DIR=getenv('CBIG_CODE_DIR'); \
addpath(fullfile(CBIG_CODE_DIR, 'external_packages', 'matlab', 'default_packages', 'cifti-matlab')); \
addpath(genpath(fullfile(CBIG_CODE_DIR, 'stable_projects', 'brain_parcellation', 'Xue2021_IndCerebellum', 'lib')));\
CBIG_IndCBM_create_surf_gifti('${surf_mesh}', '${script_dir}'); \
exit;"
fi

lh_gifti="${script_dir}/L.${surf_mesh}.func.gii"
rh_gifti="${script_dir}/R.${surf_mesh}.func.gii"
cmd="wb_command -cifti-create-dense-scalar ${template} -volume ${mask_file} ${label_mask_file}"
cmd="${cmd} -left-metric ${lh_gifti} -right-metric ${rh_gifti}"
eval $cmd
echo "Done generating template for ${name}."

rm -f ${lh_file}
rm -f ${rh_file}
rm -f ${lh_gifti}
rm -f ${rh_gifti}
rm -f ${label_mask_file}
rm -f ${mask_file}
