#!/bin/csh
# Written by Christopher Lin and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set fcount = $1
set frame_dir = $2
set MNI_ref_id = $3
set temp_2mm = $4
set LF = $5

set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`
set input = $frame_dir/${fcount_str}_MNI1mm.nii.gz
set output = $frame_dir/${fcount_str}_MNI1mm_MNI2mm.nii.gz

if(-e $output) then
    echo "[native2mni]: $output already exists." |& tee -a $LF;
    exit 0;
else
    set cmd = (mri_vol2vol --mov $input --s $MNI_ref_id --targ $temp_2mm --o $output --regheader --no-save-reg)
    echo $cmd |& tee -a $LF
    eval $cmd
    if(-e $output) then
        echo "    [native2mni]: downsample to $output finished." |& tee -a $LF
    else
        echo "    ERROR: downsample to $output failed." |& tee -a $LF
        exit 1;
    endif
endif