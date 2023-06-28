#!/bin/bash
#####
# This script preprocesses T1 and diffusion data for tractography. Parcellations and tractograms are generated. 
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####
###############
# set up environment
###############
sub=$1
algo=$2
scriptdir=$3
input_dir=$4
output_dir=$5
py_env=$6
mask_output_dir=$7
tract_streamlines=$8

# read parcellations to use
IFS=', ' read -r -a parcellation_arr <<< $9
IFS=', ' read -r -a parcels_num_arr <<< ${10}

label_dir=$scriptdir/labels
atlas_dir=$CBIG_CODE_DIR/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations
logdir=$output_dir/$algo/logs
sub_outdir=$output_dir/$algo/output/$sub
diff_dir=$input_dir/diffusion/$sub
export SUBJECTS_DIR=$input_dir/recon_all

if [ ! -d $sub_outdir ]; then mkdir -p $sub_outdir; fi
if [ ! -d $logdir ]; then mkdir -p $logdir; fi
if [ ! -d $SUBJECTS_DIR ]; then mkdir -p $SUBJECTS_DIR; fi

###############
# create symbolic links
###############
echo "[INFO]: Retrieving parcellation files"
# create symlink
echo "Step 1: Check if required files exist"
curr_dir=$( pwd )
cd $SUBJECTS_DIR
# NOTE: Please ensure annot files are in fsaverage space in order for the 
# conversion to work properly
if [ ! -d $SUBJECTS_DIR/fsaverage ]; then 
    mkdir -p fsaverage
    cd fsaverage
    ln -s $FREESURFER_HOME/subjects/fsaverage/mri
    ln -s $FREESURFER_HOME/subjects/fsaverage/mri.2mm
    ln -s $FREESURFER_HOME/subjects/fsaverage/scripts
    ln -s $FREESURFER_HOME/subjects/fsaverage/surf 
    mkdir $SUBJECTS_DIR/fsaverage/label
    cd label
    ln -s $FREESURFER_HOME/subjects/fsaverage/label/* .
    cd $curr_dir
fi
cd $curr_dir

###############
# process diffusion images
###############
echo "[INFO]: Starting diffusion image processing"
source activate $py_env
# decompress MR images
echo "Step 2: decompress MR images"
sub_dwi=$diff_dir/${sub}.nii.gz
if [ ! -e $sub_dwi ]; then
    echo "nii.gz file not found, trying nii format"
    sub_dwi=$diff_dir/${sub}.nii
    echo $sub_dwi
    if [ ! -e $sub_dwi ]; then
        echo "ERROR: Diffusion file was neither in nii.gz nor nii format"
        exit
    fi
fi
mrconvert $sub_dwi $sub_outdir/DWI.mif -fslgrad $diff_dir/${sub}.bvec $diff_dir/${sub}.bval \
    -datatype float32 -strides 0,0,0,1

# generate b0 mask
if [[ $mask_output_dir == 'NIL' ]]; then
    mask_output_dir=$output_dir/$algo/output/b0_mask
    if [ ! -d  $mask_output_dir/$sub ]; then mkdir -p $mask_output_dir/$sub; fi
    echo "Generating brain mask from b=0 image..."
    if [ ! -e  $mask_output_dir/$sub/${sub}_bet_b0_mask.nii.gz ]; then
        $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/utilities/CBIG_DiffProc_generate_b0mask.sh \
            $diff_dir $sub $mask_output_dir
    else
        echo "Brain mask already exists, skipping..."
    fi
fi
mask=$mask_output_dir/${sub}/${sub}_bet_b0_mask.nii.gz

# generate mean b0 and register to T1
echo "Step 3: generate mean b0 and perform registrations"
dwiextract $sub_outdir/DWI.mif - -bzero | mrmath - mean $sub_outdir/meanb0.mif -axis 3
mrconvert $sub_outdir/meanb0.mif $sub_outdir/meanb0.nii.gz
mrconvert $SUBJECTS_DIR/$sub/mri/brain.mgz $sub_outdir/brain.nii.gz
flirt -in $sub_outdir/meanb0.nii.gz -ref $sub_outdir/brain.nii.gz -dof 6 -omat $sub_outdir/diff2struct_fsl.mat
transformconvert $sub_outdir/diff2struct_fsl.mat $sub_outdir/meanb0.nii.gz $sub_outdir/brain.nii.gz \
    flirt_import $sub_outdir/diff2struct_mrtrix.txt

###############
# generate parcellations
###############
# generate volume file from annot
echo "Step 4: Generating parcellation from fsaverage"
len=${#parcellation_arr[@]}
for (( i=0; i<$len; i++ )); do
    parcellation=${parcellation_arr[i]}
    parcels=${parcels_num_arr[i]}
    echo " Current parcellation: [ $parcellation ]"
    echo " Number of parcels: [ $parcels ]"
    cmd="$scriptdir/CBIG_DiffProc_tractography_1_gen_parc.sh $parcellation $parcels \
    $SUBJECTS_DIR $sub $atlas_dir $sub_outdir $label_dir"
    ssh headnode "source activate $py_env; \
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd '$cmd' -walltime 02:00:00 \
    -name 'gen_parcellation' -mem 18GB -joberr '$logdir' -jobout '$logdir'" < /dev/null
done
# do label conversion for DK atlas
DK_parc='DK_82Parcels'
if [ ! -d $sub_outdir/$DK_parc ]; then mkdir -p $sub_outdir/$DK_parc; fi
labelconvert $SUBJECTS_DIR/$sub/mri/aparc+aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt \
    $label_dir/fs_DK_82.txt $sub_outdir/$DK_parc/${DK_parc}_unreg.mif
mrtransform $sub_outdir/$DK_parc/${DK_parc}_unreg.mif -linear $sub_outdir/diff2struct_mrtrix.txt \
    -inverse $sub_outdir/$DK_parc/${DK_parc}.mif
# generate tissue segmented image
echo "Step 5: generate tissue segmented image"
5ttgen freesurfer $SUBJECTS_DIR/$sub/mri/aparc+aseg.mgz $sub_outdir/5TT_unreg.mif -scratch $sub_outdir
mrtransform $sub_outdir/5TT_unreg.mif -linear $sub_outdir/diff2struct_mrtrix.txt -inverse $sub_outdir/5TT.mif

###############
# generate tractogram
###############
echo "[INFO]: Starting tractogram generation"
echo "Step 6: estimate response function and generate tractogram"
if [ ! -e $sub_outdir/dwi_wm_weights.csv ]; then
    cmd="$scriptdir/CBIG_DiffProc_tractography_2_gen_tractogram.sh \
        $sub $sub_outdir $diff_dir $algo $mask $tract_streamlines"
    ssh headnode "source activate $py_env; \
    $CBIG_CODE_DIR/setup/CBIG_pbsubmit -cmd '$cmd' -walltime 04:30:00 -ncpus 8 \
    -name 'gen_tractogram' -mem 6GB -joberr '$logdir' -jobout '$logdir'" < /dev/null
else
    echo "[WARNING]: SIFT2 output file exists! Skipping tractogram job submission..."
fi
