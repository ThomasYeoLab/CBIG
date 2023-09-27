#!/bin/bash
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

project_dir=${1} # output directory
data_type=${2}   # data_type = metric/label. 
fsLR_mesh=${3}
version=${4}     # version = 20160827 or 20170508 
fsLR_surface_dir=$CBIG_CODE_DIR/data/templates/surface

#############################################
## Setting output file extension
#############################################

if [ "$data_type" == "label" ]; then
    # To project label data, output file extension is .label.gii         
    in_extension=label.gii
elif [ "$data_type" == "metric" ]; then
    # To project metric data, output file extension is .func.gii
    in_extension=func.gii     
else
    echo "Unrecognized data_type = $data_type ! data_type = label/metric."
fi

for hemi in {L,R}; do

    if [ "$fsLR_mesh" == "fs_LR_32k" ]; then

        #############################################
        ## Upsampling from fsLR_32k to fsLR_164k
        #############################################

        echo "fsLR_mesh = $fsLR_mesh, upsample to fs_LR_164k!"

        upsample_in=$project_dir/${data_type}_${hemi}.$fsLR_mesh.$in_extension
        orig_sphere=$fsLR_surface_dir/fs_LR_32k/cifti/standard.$hemi.sphere.32k_fs_LR.surf.gii
        upsample_sphere=$fsLR_surface_dir/fs_LR_164k/cifti/standard.$hemi.sphere.164k_fs_LR.surf.gii
        upsample_out=/$project_dir/${data_type}_$hemi.fs_LR_164k.$in_extension
        orig_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.32k_fs_LR.shape.gii
        upsample_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.164k_fs_LR.shape.gii

        wb_command -${data_type}-resample $upsample_in $orig_sphere $upsample_sphere ADAP_BARY_AREA $upsample_out -area-metrics $orig_area $upsample_area

    elif [ "$fsLR_mesh" == "fs_LR_164k" ]; then

        echo "fsLR_mesh = $fsLR_mesh, no need to do upsampling!"

    else

        echo "Unrecognized fsLR_mesh = $fsLR_mesh ! fsLR_mesh = fs_LR_32k/fs_LR_164k."

    fi

    #############################################
    ## Projecting from fsLR_164k to fsaverage
    #############################################

    data_in=$project_dir/${data_type}_${hemi}.fs_LR_164k.$in_extension
    curr_sphere=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR-deformed_to-fsaverage.$hemi.sphere.164k_fs_LR.surf.gii
    new_sphere=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fsaverage_std_sphere.$hemi.164k_fsavg_$hemi.surf.gii

    data_out=$project_dir/${data_type}_${hemi}_gifti.gii
    curr_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fs_LR.$hemi.midthickness_va_avg.164k_fs_LR.shape.gii
    new_area=$fsLR_surface_dir/standard_mesh_atlases_${version}/resample_fsaverage/fsaverage.$hemi.midthickness_va_avg.164k_fsavg_$hemi.shape.gii

    wb_command -${data_type}-resample $data_in $curr_sphere $new_sphere ADAP_BARY_AREA $data_out -area-metrics $curr_area $new_area

done;

