### assumes you have the $fsav subject in your exampledir, might be changed when we use freesurfer 5
#!/bin/bash
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

workingDir=${1}
segment_name=${2}
fsav=fsaverage

export SUBJECTS_DIR=${workingDir}

# prepare subject
subjectDir=${workingDir}/${subName}
echo $subjectDir
left_map=${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/data/fs_L/
right_map=${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/data/fs_R/


# Convert to Fsaverage
cd $right_map
echo right hemisphere
caret_command -deformation-map-apply ${right_map}/fsaverage.LR.registered-to-fs_R_fix.164k_fs_R.deform_map METRIC_NEAREST_NODE ${subjectDir}/${segment_name}.R.164k_fs_LR.label.gii ${subjectDir}/${segment_name}.R.fs_average.func.gii

cd $left_map
echo left hemisphere
caret_command -deformation-map-apply ${left_map}/fsaverage.LR.registered-to-fs_L_fix.164k_fs_L.deform_map METRIC_NEAREST_NODE ${subjectDir}/${segment_name}.L.164k_fs_LR.label.gii ${subjectDir}/${segment_name}.L.fs_average.func.gii

