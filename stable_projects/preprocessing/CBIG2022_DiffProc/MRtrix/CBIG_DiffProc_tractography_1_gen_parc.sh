#!/bin/bash
#####
# This script projects parcellations to the native volume space and converts them to a mif format. 
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####
###############
# set up environment
###############
parcellation=$1
parcels=$2
export SUBJECTS_DIR=$3
sub=$4
atlas_dir=$5
sub_outdir=$6
label_dir=$7

###############
# generate volume and mif for parcellation
###############
echo " Current parcellation: [ $parcellation ]"
echo " Number of parcels: [ $parcels ]"
# project to subject space
lh_fsaverage_annot=$atlas_dir/FreeSurfer5.3/fsaverage/label/lh.${parcellation}_order.annot
rh_fsaverage_annot=$atlas_dir/FreeSurfer5.3/fsaverage/label/rh.${parcellation}_order.annot

if [ ! -e $SUBJECTS_DIR/$sub/label/lh.aparc.$parcellation.annot ]; then
	echo "Projecting parcellation to individual T1 space (lh)"
	mri_surf2surf --srcsubject fsaverage --trgsubject $sub --hemi lh \
	--sval-annot $lh_fsaverage_annot \
	--tval $SUBJECTS_DIR/$sub/label/lh.aparc.$parcellation.annot
fi
if [ ! -e $SUBJECTS_DIR/$sub/label/rh.aparc.$parcellation.annot ]; then
	echo "Projecting parcellation to individual T1 space (rh)"
	mri_surf2surf --srcsubject fsaverage --trgsubject $sub --hemi rh \
	--sval-annot $rh_fsaverage_annot \
	--tval $SUBJECTS_DIR/$sub/label/rh.aparc.$parcellation.annot 
fi
# generate volume (mgz)
if [ ! -e $sub_outdir/$parcellation ]; then mkdir -p $sub_outdir/$parcellation; fi
if [ ! -e $sub_outdir/$parcellation/${parcellation}_lh.mgz ]; then
	mri_label2vol --annot $SUBJECTS_DIR/$sub/label/lh.aparc.$parcellation.annot --subject $sub --hemi lh \
		--temp $SUBJECTS_DIR/$sub/mri/aparc+aseg.mgz --identity --o $sub_outdir/$parcellation/${parcellation}_lh.mgz
fi
if [ ! -e $sub_outdir/$parcellation/${parcellation}_rh.mgz ]; then
	mri_label2vol --annot $SUBJECTS_DIR/$sub/label/rh.aparc.$parcellation.annot --subject $sub --hemi rh \
		--temp $SUBJECTS_DIR/$sub/mri/aparc+aseg.mgz --identity --o $sub_outdir/$parcellation/${parcellation}_rh.mgz
fi
# combine labels from both hemispheres
mrconvert $sub_outdir/$parcellation/${parcellation}_lh.mgz $sub_outdir/$parcellation/${parcellation}_lh.nii.gz -force
mrconvert $sub_outdir/$parcellation/${parcellation}_rh.mgz $sub_outdir/$parcellation/${parcellation}_rh.nii.gz -force
fslmaths $sub_outdir/$parcellation/${parcellation}_lh.nii.gz -bin -sub 1 -mul -1 \
	$sub_outdir/$parcellation/${parcellation}_lh_mod.nii.gz 
# add half the parcel counts to rh volume
fslmaths $sub_outdir/$parcellation/${parcellation}_rh.nii.gz \
	-mul $sub_outdir/$parcellation/${parcellation}_lh_mod.nii.gz -add $(($parcels / 2)) \
	-thr $(($parcels / 2 + 1)) $sub_outdir/$parcellation/${parcellation}_rh_mod.nii.gz -odt int
mri_concat --combine --i $sub_outdir/$parcellation/${parcellation}_lh.mgz \
	--i $sub_outdir/$parcellation/${parcellation}_rh_mod.nii.gz \
	--o $sub_outdir/$parcellation/${parcellation}_combined_nosubcort.mgz
# add subcortical volumes
if [ ! -e $sub_outdir/$parcellation/fs_subcort.nii.gz ]; then
labelconvert $SUBJECTS_DIR/$sub/mri/aparc+aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt \
	$label_dir/fs_subcort.txt $sub_outdir/$parcellation/fs_subcort.mif
mrconvert $sub_outdir/$parcellation/fs_subcort.mif $sub_outdir/$parcellation/fs_subcort.nii.gz -force
fi
fslmaths $sub_outdir/$parcellation/fs_subcort.nii.gz -add $parcels \
	-thr $(($parcels + 1)) $sub_outdir/$parcellation/${parcellation}_subcort_mod.nii.gz \
	-odt int # add the parcel counts to subcort volume
mri_concat --combine --i $sub_outdir/$parcellation/${parcellation}_combined_nosubcort.mgz \
	--i $sub_outdir/$parcellation/${parcellation}_subcort_mod.nii.gz \
	--o $sub_outdir/$parcellation/${parcellation}_combined.mgz

# convert to mrtrix format
echo "Converting parcellation to MRtrix format"
mrconvert -datatype Uint32 -force $sub_outdir/$parcellation/${parcellation}_combined.mgz \
	$sub_outdir/$parcellation/${parcellation}_unreg.mif
mrtransform $sub_outdir/$parcellation/${parcellation}_unreg.mif -linear $sub_outdir/diff2struct_mrtrix.txt \
	-inverse $sub_outdir/$parcellation/${parcellation}.mif
