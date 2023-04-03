#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Yan2023_homotopic
# remove useless stable projects
rm -r Standalone_Yan2023_homotopic/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering
rm -r Standalone_Yan2023_homotopic/stable_projects/brain_parcellation/Xue2021_IndCerebellum
rm -r Standalone_Yan2023_homotopic/stable_projects/brain_parcellation/Kong2022_ArealMSHBM
rm -r Standalone_Yan2023_homotopic/stable_projects/brain_parcellation/Kong2019_MSHBM
rm -r Standalone_Yan2023_homotopic/stable_projects/disorder_subtypes
rm -r Standalone_Yan2023_homotopic/stable_projects/fMRI_dynamics
rm -r Standalone_Yan2023_homotopic/stable_projects/meta-analysis
rm -r Standalone_Yan2023_homotopic/stable_projects/predict_phenotypes
rm -r Standalone_Yan2023_homotopic/stable_projects/preprocessing
rm -r Standalone_Yan2023_homotopic/stable_projects/registration
