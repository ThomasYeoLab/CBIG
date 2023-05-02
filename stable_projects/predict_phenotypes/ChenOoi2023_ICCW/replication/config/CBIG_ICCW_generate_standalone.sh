#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_ChenOoi2023_ICCW
# remove unused stable projects
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/brain_parcellation
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/disorder_subtypes
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/fMRI_dynamics
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/registration
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/meta-analysis
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/preprocessing
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/predict_phenotypes/He2019_KRDNN
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/predict_phenotypes/Nguyen2020_RNNAD
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/predict_phenotypes/Ooi2022_MMP
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/predict_phenotypes/An2022_gcVAE
rm -r Standalone_ChenOoi2023_ICCW/stable_projects/predict_phenotypes/Kong2023_GradPar