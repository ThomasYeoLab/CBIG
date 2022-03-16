#!/bin/sh
#
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git $CBIG_CODE_DIR/* Standalone_He2022_MM

# remove useless stable projects
rm -r Standalone_He2022_MM/stable_projects/brain_parcellation
rm -r Standalone_He2022_MM/stable_projects/disorder_subtypes
rm -r Standalone_He2022_MM/stable_projects/fMRI_dynamics
rm -r Standalone_He2022_MM/stable_projects/meta-analysis
rm -r Standalone_He2022_MM/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
rm -r Standalone_He2022_MM/stable_projects/predict_phenotypes/He2019_KRDNN
rm -r Standalone_He2022_MM/stable_projects/predict_phenotypes/Nguyen2020_RNNAD
rm -r Standalone_He2022_MM/stable_projects/predict_phenotypes/Ooi2022_MMP
rm -r Standalone_He2022_MM/stable_projects/preprocessing/CBIG_fMRI_Preproc2016
rm -r Standalone_He2022_MM/stable_projects/preprocessing/CBIG2022_DiffProc
rm -r Standalone_He2022_MM/stable_projects/preprocessing/Orban2020_tod
rm -r Standalone_He2022_MM/stable_projects/registration