#!/bin/sh
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Kong2021_pMFM
# remove useless stable projects
rm -r Standalone_Kong2021_pMFM/stable_projects/fMRI_dynamics/Liegeois2017_Surrogates
rm -r Standalone_Kong2021_pMFM/stable_projects/fMRI_dynamics/Wang2018_MFMem
rm -r Standalone_Kong2021_pMFM/stable_projects/brain_parcellation
rm -r Standalone_Kong2021_pMFM/stable_projects/disorder_subtypes
rm -r Standalone_Kong2021_pMFM/stable_projects/registration
rm -r Standalone_Kong2021_pMFM/stable_projects/meta-analysis
rm -r Standalone_Kong2021_pMFM/stable_projects/preprocessing
rm -r Standalone_Kong2021_pMFM/stable_projects/predict_phenotypes