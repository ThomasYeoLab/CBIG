### assumes you have the $fsav subject in your exampledir, might be changed when we use freesurfer 5
#!/bin/bash
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

#wb=/data/users/aschaefer/src/workbench_centos/workbench/bin_rh_linux64/wb_command
workingDir=${1}
segment_name=${2}
type_of_smoothing=${3}
fsav=fsaverage
inter_atlas_deformation_fslr=${CBIG_CODE_DIR}/utilities/matlab/fslr_matlab/data/deform_map/
freesurferDir=$workingDir ### working freesurfer Dir
export SUBJECTS_DIR=${freesurferDir}

# prepare subject
subjectDir=${workingDir}
echo $subjectDir



#### aplly registration
#mri_surf2surf --srcsubject est --srcsurfval ${subjectDir}/lh_borders.mgh --trgsurfval ${subjectDir}/test_lh.mgh --hemi lh --trgsubject fsaverage --noreshape
mri_convert ${subjectDir}/${segment_name}_lh_borders.mgh ${subjectDir}/${segment_name}_test_lh.gii
mri_convert ${subjectDir}/${segment_name}_rh_borders.mgh ${subjectDir}/${segment_name}_test_rh.gii
echo mri_convert ${subjectDir}/${segment_name}_lh_borders.mgh ${subjectDir}/${segment_name}_test_lh.gii

# Convert to FS_LR

cd ${inter_atlas_deformation_fslr}
ls

caret_command -deformation-map-apply fsaverage.R.registered-to-fs_LR_fix.164k_fs_LR.deform_map $type_of_smoothing ${subjectDir}/${segment_name}_test_rh.gii  ${subjectDir}/${segment_name}.R.fs_LR.func.gii

caret_command -deformation-map-apply fsaverage.L.registered-to-fs_LR_fix.164k_fs_LR.deform_map ${type_of_smoothing} ${subjectDir}/${segment_name}_test_lh.gii ${subjectDir}/${segment_name}.L.fs_LR.func.gii

cd ${subjectDir}/${subName}/
