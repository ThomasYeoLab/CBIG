#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_ChenTam2022_TRBPC
# remove useless stable projects
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/brain_parcellation
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/fMRI_dynamics
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/registration
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/meta-analysis
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/preprocessing
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/predict_phenotypes/He2019_KRDNN
rm -r Standalone_ChenTam2022_TRBPC/stable_projects/predict_phenotypes/Nguyen2020_RNNAD
