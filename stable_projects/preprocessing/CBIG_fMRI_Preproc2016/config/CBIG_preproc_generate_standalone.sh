#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_CBIG_fMRI_Preproc2016
# remove useless stable projects
rm -r Standalone_CBIG_fMRI_Preproc2016/stable_projects/brain_parcellation/Kong2019_MSHBM
rm -r Standalone_CBIG_fMRI_Preproc2016/stable_projects/disorder_subtypes
rm -r Standalone_CBIG_fMRI_Preproc2016/stable_projects/fMRI_dynamics
rm -r Standalone_CBIG_fMRI_Preproc2016/stable_projects/registration
