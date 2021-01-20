#!/bin/bash
# Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Sun2019_ADJointFactors
# remove useless stable projects
rm -r Standalone_Sun2019_ADJointFactors/stable_projects/brain_parcellation
rm -r Standalone_Sun2019_ADJointFactors/stable_projects/disorder_subtypes/Tang2019_ASDFactors
rm -r Standalone_Sun2019_ADJointFactors/stable_projects/fMRI_dynamics
rm -r Standalone_Sun2019_ADJointFactors/stable_projects/preprocessing
rm -r Standalone_Sun2019_ADJointFactors/stable_projects/registration
