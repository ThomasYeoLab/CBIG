#!/bin/bash
# Written by Siyi Tang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Tang2020_ASDFactors
# remove useless stable projects
rm -r Standalone_Tang2020_ASDFactors/stable_projects/brain_parcellation/Kong2019_MSHBM
rm -r Standalone_Tang2020_ASDFactors/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering
rm -r Standalone_Tang2020_ASDFactors/stable_projects/disorder_subtypes/Kebets2019_TransdiagnosticComponents
rm -r Standalone_Tang2020_ASDFactors/stable_projects/disorder_subtypes/Sun2019_ADJointFactors
rm -r Standalone_Tang2020_ASDFactors/stable_projects/fMRI_dynamics
rm -r Standalone_Tang2020_ASDFactors/stable_projects/meta-analysis
rm -r Standalone_Tang2020_ASDFactors/stable_projects/predict_phenotypes
rm -r Standalone_Tang2020_ASDFactors/stable_projects/preprocessing
rm -r Standalone_Tang2020_ASDFactors/stable_projects/registration
