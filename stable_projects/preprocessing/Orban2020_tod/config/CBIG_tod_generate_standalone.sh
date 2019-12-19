#!/bin/bash
# Written by Csaba Orban and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Orban2020_time_of_day
# remove useless stable projects
rm -r Standalone_Orban2020_time_of_day/stable_projects/brain_parcellation
rm -r Standalone_Orban2020_time_of_day/stable_projects/disorder_subtypes
rm -r Standalone_Orban2020_time_of_day/stable_projects/fMRI_dynamics
rm -r Standalone_Orban2020_time_of_day/stable_projects/preprocessing/CBIG_fMRI_Preproc2016 
rm -r Standalone_Orban2020_time_of_day/stable_projects/preprocessing/Li2019_GSR
rm -r Standalone_Orban2020_time_of_day/stable_projects/registration
rm -r Standalone_Orban2020_time_of_day/stable_projects/meta-analysis
rm -r Standalone_Orban2020_time_of_day/stable_projects/predict_phenotypes
