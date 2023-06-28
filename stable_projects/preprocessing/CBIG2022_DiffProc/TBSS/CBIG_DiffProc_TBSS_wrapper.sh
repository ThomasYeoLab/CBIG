#!/bin/bash

#####
# Example: 
#    $CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/ \
#        TBSS/CBIG_DiffProc_TBSS_wrapper.sh path/to/tbss/dir
#
# This function runs the TBSS workflow, assuming the file directory has already been set up. 
# Refer to the readme in the TBSS folder for instructions.
#
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
#####

### Export directories
export scriptdir=$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG2022_DiffProc/TBSS
export tbssdir=$1
export FAdir=$tbssdir/FA

### Run TBSS scripts with default settings
cd $tbssdir
echo "------------ Run TBSS workflow ------------"
echo "[1] tbss_1_preproc: Move FA images to single folder and create masks"
$scriptdir/tbss_1_preproc_torque *.nii.gz
echo "[2] tbss_2_reg: Register subject FA images to template"
$scriptdir/tbss_2_reg_torque -T
echo "[3] tbss_3_postreg: Apply registrations and create FA skeleton"
$scriptdir/tbss_3_postreg_torque -S
echo "[4] tbss_4_prestats: Threshold FA skeleton"
tbss_4_prestats 0.2
### optional step
echo "[5] tbss_non_FA: Apply saved transformation from FA to other metrics and calculate skeletons" 
$scriptdir/tbss_non_FA_torque MD