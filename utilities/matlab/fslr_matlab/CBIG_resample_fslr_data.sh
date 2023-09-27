#!/bin/bash
# Written by Ruby Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

orig_mesh=${1} # original data mesh
resample_mesh=${2} # target mesh of resampling
project_dir=${3} # output directory
data_type=${4}   # data_type = metric/label.
version=${5}     # version = 20160827 or 20170508 
fsLR_surface_dir=$CBIG_CODE_DIR/data/templates/surface


#############################################
## Resampling between fsLR_164k and fsLR_32k
#############################################
if [ "$orig_mesh" == "fs_LR_32k" ]; then
    standard_orig_mesh=32k_fs_LR;
fi
if [ "$orig_mesh" == "fs_LR_164k" ]; then
    standard_orig_mesh=164k_fs_LR;
fi
if [ "$resample_mesh" == "fs_LR_32k" ]; then
    standard_resample_mesh=32k_fs_LR;
fi
if [ "$resample_mesh" == "fs_LR_164k" ]; then
    standard_resample_mesh=164k_fs_LR;
fi

if [ "$data_type" == "metric" ]; then
    fsLR_extension=func;
fi
if [ "$data_type" == "label" ]; then
    fsLR_extension=label;
fi

for hemi in {L,R}; do
    resample_in=${project_dir}/${data_type}_${hemi}.${orig_mesh}.${fsLR_extension}.gii
    resample_sphere=$fsLR_surface_dir/${resample_mesh}/cifti/standard.$hemi.sphere.${standard_resample_mesh}.surf.gii
    orig_sphere=$fsLR_surface_dir/${orig_mesh}/cifti/standard.$hemi.sphere.${standard_orig_mesh}.surf.gii
    resample_out=${project_dir}/${data_type}_${hemi}.${resample_mesh}.${fsLR_extension}.gii
    resample_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.${standard_resample_mesh}.shape.gii
    orig_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.${standard_orig_mesh}.shape.gii

    wb_command -${data_type}-resample $resample_in $orig_sphere $resample_sphere ADAP_BARY_AREA $resample_out -area-metrics $orig_area $resample_area
done
