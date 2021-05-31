#!/bin/sh
# Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Kong2022_ArealMSHBM
# remove useless stable projects
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/disorder_subtypes
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/fMRI_dynamics
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/registration
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/meta-analysis
rm -r Standalone_Kong2022_ArealMSHBM/stable_projects/preprocessing
