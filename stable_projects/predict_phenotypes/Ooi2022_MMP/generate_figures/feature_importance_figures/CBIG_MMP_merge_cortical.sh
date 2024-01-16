#!/bin/bash
#####
# This script combines images generated from freesurfer.
# 1. lh_lateral: lateral freesurfer image for left hemisphere
# 2. rh_lateral: lateral freesurfer image for right hemisphere
# 3. lh_medial:  medial freesurfer image for left hemisphere
# 4. rh_medial: medial freesurfer image for right hemisphere
# 5. tmp_dir: temporary directory where files can be combined
# 6. output: output name for the final file
# 
# EXAMPLE: 
#    CBIG_MMP_merge_cortical.sh $lh_lateral $rh_lateral $lh_medial $rh_medial \
#        $tmp_dir #output
# Written by Leon Ooi and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# read in file names
lh_lateral=$1
rh_lateral=$2
lh_medial=$3
rh_medial=$4
tmp_dir=$5
output=$tmp_dir/$6

# combine hemispheres
convert -background white $lh_lateral $rh_lateral +append $tmp_dir/tmp_lateral.png
convert -background white $lh_medial $rh_medial +append $tmp_dir/tmp_medial.png
# combine sides
convert -background white $tmp_dir/tmp_lateral.png $tmp_dir/tmp_medial.png -append $output
# remove temporary files
rm $tmp_dir/tmp_lateral.png $tmp_dir/tmp_medial.png
