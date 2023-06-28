#!/bin/bash
#####
# This script generates the connectome from the tractogram. Additionally fits a DTI model and generates the 
# average of DTI indices over each connection from ROI A to ROI B.
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####
###############
# set up environment
###############
sub=$1
algo=$2
input_dir=$3
output_dir=$4
tract_streamlines=$5
sub_outdir=$output_dir/$algo/output/$sub
export SUBJECTS_DIR=$input_dir/recon_all
if [ ! -d $sub_outdir/connectomes ]; then mkdir -p $sub_outdir/connectomes; fi

# read parcellations to use
IFS=', ' read -r -a parcellation_arr <<< $6
IFS=', ' read -r -a parcels_num_arr <<< $7

###############
# generate connectome
###############
for parcellation in ${parcellation_arr[@]}; do
    echo "[PROGRESS]: Current connectome = [ $sub_outdir/connectomes/connectome_${parcellation}_SIFT2.csv ]"
    tck2connectome -tck_weights_in $sub_outdir/dwi_wm_weights.csv -assignment_radial_search 2 \
        -stat_edge sum -symmetric $sub_outdir/${tract_streamlines}.tck \
        $sub_outdir/${parcellation}/${parcellation}.mif \
        $sub_outdir/connectomes/connectome_${parcellation}_SIFT2.csv 
    echo "[PROGRESS]: Current connectome = [ $sub_outdir/connectomes/connectome_${parcellation}_unfiltered.csv ]"
    tck2connectome -assignment_radial_search 2 -stat_edge sum -symmetric \
        $sub_outdir/${tract_streamlines}.tck $sub_outdir/${parcellation}/${parcellation}.mif \
        $sub_outdir/connectomes/connectome_${parcellation}_unfiltered.csv -force
    echo "[PROGRESS]: Current connectome = [ $sub_outdir/connectomes/connectome_${parcellation}_length.csv ]"
    tck2connectome -assignment_radial_search 2 -stat_edge mean -symmetric -scale_length \
        $sub_outdir/${tract_streamlines}.tck $sub_outdir/${parcellation}/${parcellation}.mif \
        $sub_outdir/connectomes/connectome_${parcellation}_length.csv 
done

###############
# extract DTI indices
###############
dwi2tensor $sub_outdir/DWI.mif $sub_outdir/DTI_tensor.mif
tensor2metric $sub_outdir/DTI_tensor.mif -fa $sub_outdir/FA.mif -adc $sub_outdir/MD.mif \
    -ad $sub_outdir/AD.mif -rd $sub_outdir/RD.mif
dti_arr=("FA" "MD" "AD" "RD")
for metric in ${dti_arr[@]}; do
    tcksample -stat_tck mean $sub_outdir/${tract_streamlines}.tck $sub_outdir/$metric.mif $sub_outdir/$metric.csv
    if [ ! -d $sub_outdir/connectomes/$metric ]; then mkdir -p $sub_outdir/connectomes/$metric; fi
    for parcellation in ${parcellation_arr[@]}; do
        echo "[PROGRESS]: Current connectome = [ $sub_outdir/connectomes/$metric/connectome_${parcellation}_$metric.csv ]"
        tck2connectome -scale_file $sub_outdir/$metric.csv -stat_edge mean -symmetric \
        $sub_outdir/${tract_streamlines}.tck $sub_outdir/${parcellation}/${parcellation}.mif \
        $sub_outdir/connectomes/$metric/connectome_${parcellation}_$metric.csv 
    done
done
