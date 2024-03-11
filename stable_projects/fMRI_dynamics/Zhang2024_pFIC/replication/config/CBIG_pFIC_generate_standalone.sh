#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Zhang2024_pFIC
# remove unused stable projects
rm -r Standalone_Zhang2024_pFIC/stable_projects/brain_parcellation
rm -r Standalone_Zhang2024_pFIC/stable_projects/disorder_subtypes
rm -r Standalone_Zhang2024_pFIC/stable_projects/fMRI_dynamics/Kong2021_pMFM
rm -r Standalone_Zhang2024_pFIC/stable_projects/fMRI_dynamics/Liegeois2017_Surrogates
rm -r Standalone_Zhang2024_pFIC/stable_projects/fMRI_dynamics/Wang2018_MFMem
rm -r Standalone_Zhang2024_pFIC/stable_projects/meta-analysis
rm -r Standalone_Zhang2024_pFIC/stable_projects/predict_phenotypes
rm -r Standalone_Zhang2024_pFIC/stable_projects/preprocessing
rm -r Standalone_Zhang2024_pFIC/stable_projects/registration

