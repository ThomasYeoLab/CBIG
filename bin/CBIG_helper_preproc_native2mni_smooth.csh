#!/bin/csh
# Written by Christopher Lin and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

set fcount = $1
set frame_dir = $2
set sm_mask = $3
set inverted_sm_mask = $4
set sm = $5
set LF = $6

set fcount_str = `echo $fcount | awk '{printf ("%04d",$1)}'`
set input = $frame_dir/${fcount_str}_MNI1mm_MNI2mm.nii.gz
mkdir -p $frame_dir/sm
set output = $frame_dir/sm/${fcount_str}_MNI1mm_MNI2mm_sm${sm}.nii.gz
#Note that fwhm = 2.35482 * std, fslmaths -s is in unit of mm, not voxel.
set std = `awk "BEGIN {print ${sm}/2.35482}"`    
if(-e $output) then
    echo "[native2mni]: $output already exists." |& tee -a $LF;
    exit 0;
else
    if($?sm_mask) then
        # if the user passes in a volume sm_mask, the procedure is
        # 1. smooth volume data within sm_mask
        # 2. smooth sm_mask within sm_mask
        # 3. divide smoothed volume by smoothed sm_mask (deal with boundary problem)
        set tmp1 = $frame_dir/tmp1_${fcount_str}.nii.gz
        set tmp2 = $frame_dir/tmp2_${fcount_str}.nii.gz
        set tmp3 = $frame_dir/tmp3_${fcount_str}.nii.gz		
        set tmp4 = $frame_dir/tmp4_${fcount_str}.nii.gz
        set input_masksmoothed = $frame_dir/input_masksmoothed_${fcount_str}.nii.gz
        set input_outsidemask_smoothed = $frame_dir/input_outsidemask_smoothed_${fcount_str}.nii.gz

        set cmd = (fslmaths $input -mas $sm_mask -s $std -mas $sm_mask $tmp1)
        
        echo $cmd |& tee -a $LF
        eval $cmd
        
        set cmd = (fslmaths $sm_mask -s $std -mas $sm_mask $tmp2)
        echo $cmd |& tee -a $LF
        eval $cmd

        set cmd = (fslmaths $tmp1 -div $tmp2 $input_masksmoothed)
            echo $cmd |& tee -a $LF
        eval $cmd

                    #4. smooth outside of the mask, but inside the brain			    
        set cmd = (fslmaths $input -mas $inverted_sm_mask -s $std -mas $inverted_sm_mask $tmp3)
        echo $cmd |& tee -a $LF
        eval $cmd

        set cmd = (fslmaths $inverted_sm_mask -s $std -mas $inverted_sm_mask $tmp4)
        echo $cmd |& tee -a $LF
        eval $cmd

        set cmd = (fslmaths $tmp3 -div $tmp4 $input_outsidemask_smoothed)
        echo $cmd |& tee -a $LF
        eval $cmd

                    #5. combine externally smoothed data with rest of the data
        set cmd = (fslmaths $input_masksmoothed -add $input_outsidemask_smoothed $output)
        echo $cmd |& tee -a $LF
        eval $cmd
    else
        set cmd = (fslmaths $input -s $std $output)
        echo $cmd |& tee -a $LF
        eval $cmd
    endif
    
    if(-e $output) then
        echo "[native2mni]: smooth to $output finished." |& tee -a $LF
    else
        echo "ERROR: smooth to $output failed." |& tee -a $LF
        exit 1;
    endif
endif