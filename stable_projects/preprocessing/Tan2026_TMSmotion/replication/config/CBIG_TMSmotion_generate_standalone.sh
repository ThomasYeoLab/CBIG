#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Tan2026_TMSmotion
# remove unrelated stable projects
rm -r Standalone_Tan2026_TMSmotion/stable_projects/disorder_subtypes
rm -r Standalone_Tan2026_TMSmotion/stable_projects/fMRI_dynamics
rm -r Standalone_Tan2026_TMSmotion/stable_projects/harmonization
rm -r Standalone_Tan2026_TMSmotion/stable_projects/meta-analysis
rm -r Standalone_Tan2026_TMSmotion/stable_projects/predict_phenotypes
rm -r Standalone_Tan2026_TMSmotion/stable_projects/registration
rm -r Standalone_Tan2026_TMSmotion/stable_projects/preprocessing/CBIG2022_DiffProc
rm -r Standalone_Tan2026_TMSmotion/stable_projects/preprocessing/Li2019_GSR
rm -r Standalone_Tan2026_TMSmotion/stable_projects/preprocessing/Orban2020_tod
rm -r Standalone_Tan2026_TMSmotion/stable_projects/brain_parcellation/Lim2026_MSHBM_epilepsy
rm -r Standalone_Tan2026_TMSmotion/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal
rm -r Standalone_Tan2026_TMSmotion/stable_projects/brain_parcellation/Xue2021_IndCerebellum
rm -r Standalone_Tan2026_TMSmotion/stable_projects/brain_parcellation/Yan2023_homotopic
rm -r Standalone_Tan2026_TMSmotion/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering