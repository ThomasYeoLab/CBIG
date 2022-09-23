#!/bin/bash
#####
# This script estimates the fibre orientation directions and generates the tractogram.
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####
###############
# set up environment
###############
sub=$1
sub_outdir=$2 
diff_dir=$3
algo=$4
mask=$5
tract_streamlines=$6

###############
# generate tractogram
###############
# estimate FODs
dwi2response msmt_5tt $sub_outdir/DWI.mif $sub_outdir/5TT.mif $sub_outdir/RF_WM.txt $sub_outdir/RF_GM.txt \
	$sub_outdir/RF_CSF.txt -scratch $sub_outdir -voxels $sub_outdir/RF_voxels.mif 
dwi2fod msmt_csd $sub_outdir/DWI.mif $sub_outdir/RF_WM.txt $sub_outdir/WM_FODs.mif \
	$sub_outdir/RF_GM.txt $sub_outdir/GM.mif $sub_outdir/RF_CSF.txt $sub_outdir/CSF.mif \
	-mask $mask

# generate tractogram
if [ $algo == 'iFOD2' ]; then
tckgen -act $sub_outdir/5TT.mif -algorithm $algo -backtrack -crop_at_gmwmi -samples 4 -nthreads 8 \
	-output_seeds $sub_outdir/out_seeds.nii.gz -power 0.330000 -seed_dynamic $sub_outdir/WM_FODs.mif \
	-select ${tract_streamlines} $sub_outdir/WM_FODs.mif $sub_outdir/${tract_streamlines}.tck
else 
tckgen -act $sub_outdir/5TT.mif -algorithm $algo -crop_at_gmwmi -samples 4 -nthreads 8 \
	-output_seeds $sub_outdir/out_seeds.nii.gz -power 0.330000 -seed_dynamic $sub_outdir/WM_FODs.mif \
	-select ${tract_streamlines} $sub_outdir/WM_FODs.mif $sub_outdir/${tract_streamlines}.tck
fi
# run SIFT2
tcksift2 -act $sub_outdir/5TT.mif -nthreads 8 \
	-out_mu $sub_outdir/dwi_wm_mu.txt $sub_outdir/${tract_streamlines}.tck \
	$sub_outdir/WM_FODs.mif $sub_outdir/dwi_wm_weights.csv
