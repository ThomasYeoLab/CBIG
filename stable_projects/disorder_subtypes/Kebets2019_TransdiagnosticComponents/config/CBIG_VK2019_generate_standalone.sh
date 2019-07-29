#!/bin/bash
# Written by Valeria Kebets and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Kebets2019_TransdiagnosticComponents
# remove useless stable projects
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/brain_parcellation
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/disorder_subtypes/Zhang2016_ADFactors
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/disorder_subtypes/Tang2019_ASDFactors
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/fMRI_dynamics
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/preprocessing
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/registration
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/meta-analysis
rm -r Standalone_Kebets2019_TransdiagnosticComponents/stable_projects/predict_phenotypes
