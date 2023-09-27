#!/bin/bash
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

project_dir=${1} # output directory
data_type=${2}   # data_type = metric/label.
version=${3}     # version = 20160827 or 20170508 
fsLR_surface_dir=$CBIG_CODE_DIR/data/templates/surface

#############################################
## Setting output file extension
#############################################

if [ "$data_type" == "label" ]; then
    # To project label data, output file extension is .label.gii         
    out_extension=label.gii
elif [ "$data_type" == "metric" ]; then
    # To project metric data, output file extension is .func.gii
    out_extension=func.gii     
else
    echo "Unrecognized data_type = $data_type ! data_type = label/metric."
fi

for hemi in {L,R}; do

    #############################################
    ## Converting from nifti format (fsaverage) to gifti format
    #############################################

    if [ ! -f $project_dir/${data_type}_${hemi}_gifti.gii ]; then
        mri_convert $project_dir/${data_type}_${hemi}_fsaverage_borders.mgh $project_dir/${data_type}_${hemi}_gifti.gii
    fi

    #############################################
    ## Projecting gifti file to fsLR_164k 
    #############################################

    data_in=$project_dir/${data_type}_${hemi}_gifti.gii
    curr_sphere=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fsaverage_std_sphere.$hemi.164k_fsavg_$hemi.surf.gii
    new_sphere=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR-deformed_to-fsaverage.$hemi.sphere.164k_fs_LR.surf.gii
    data_out=/$project_dir/${data_type}_$hemi.fs_LR_164k.$out_extension
    curr_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fsaverage.$hemi.midthickness_va_avg.164k_fsavg_$hemi.shape.gii
    new_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.164k_fs_LR.shape.gii

    wb_command -${data_type}-resample $data_in $curr_sphere $new_sphere ADAP_BARY_AREA $data_out -area-metrics $curr_area $new_area

    #############################################
    ## Downsampling from fsLR_164k to fsLR_32k
    #############################################

    downsample_in=$data_out
    orig_sphere=$fsLR_surface_dir/fs_LR_164k/cifti/standard.$hemi.sphere.164k_fs_LR.surf.gii
    downsample_sphere=$fsLR_surface_dir/fs_LR_32k/cifti/standard.$hemi.sphere.32k_fs_LR.surf.gii
    downsample_out=/$project_dir/${data_type}_$hemi.fs_LR_32k.$out_extension
    orig_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.164k_fs_LR.shape.gii
    downsample_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.32k_fs_LR.shape.gii

    wb_command -${data_type}-resample $downsample_in $orig_sphere $downsample_sphere ADAP_BARY_AREA $downsample_out -area-metrics $orig_area $downsample_area

done;

