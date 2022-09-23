#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_CBIG2022_DiffProc

# remove useless stable projects
rm -r Standalone_CBIG2022_DiffProc/stable_projects/brain_parcellation
rm -r Standalone_CBIG2022_DiffProc/stable_projects/disorder_subtypes
rm -r Standalone_CBIG2022_DiffProc/stable_projects/fMRI_dynamics
rm -r Standalone_CBIG2022_DiffProc/stable_projects/registration
rm -r Standalone_CBIG2022_DiffProc/stable_projects/predict_phenotypes
rm -r Standalone_CBIG2022_DiffProc/stable_projects/meta-analysis
rm -r Standalone_CBIG2022_DiffProc/stable_projects/preprocessing/Orban2020_tod
rm -r Standalone_CBIG2022_DiffProc/stable_projects/preprocessing/Li2019_GSR
rm -r Standalone_CBIG2022_DiffProc/stable_projects/preprocessing/CBIG_fMRI_Preproc2016