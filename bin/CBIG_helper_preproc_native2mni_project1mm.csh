#!/bin/csh

# Written by Christopher Lin and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
set fcount = $1
set frame_dir = $2
# set input $3
# set output $4

set anat_s = $3
set anat_dir = $4
set MNI_ref_id = $5
set MNI_ref_dir = $6
set regfile = $7
set LF = $8

echo FRAME NUMBER $fcount |& tee -a $LF

set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`
		
set input = $frame_dir/orig_frames${fcount_str}.nii.gz
set output = $frame_dir/${fcount_str}_MNI1mm.nii.gz
if(-e $output) then
    echo "    [native2mni]: $output already exists." |& tee -a $LF;
    exit 0;
else
    set cmd = (CBIG_vol2vol_m3z.csh -src-id $anat_s -src-dir $anat_dir -targ-id $MNI_ref_id \
    -targ-dir $MNI_ref_dir -in $input -out $output -reg $regfile -no-cleanup)
    echo $cmd |& tee -a $LF
    eval $cmd |& tee -a $LF
    if(-e $output) then
        echo "    [native2mni]: projection to $output finished." |& tee -a $LF
    else
        echo "    ERROR: projection to $output failed." |& tee -a $LF
        exit 1;
    endif
endif