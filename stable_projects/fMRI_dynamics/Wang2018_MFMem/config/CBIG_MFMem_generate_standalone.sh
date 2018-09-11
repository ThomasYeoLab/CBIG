#!/bin/bash
# Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git CBIG/* Standalone_Wang2018_MFMem
# remove useless stable projects
rm -r Standalone_Wang2018_MFMem/stable_projects/brain_parcellation
rm -r Standalone_Wang2018_MFMem/stable_projects/disorder_subtypes
rm -r Standalone_Wang2018_MFMem/stable_projects/fMRI_dynamics/Liegeois2017_Surrogates
rm -r Standalone_Wang2018_MFMem/stable_projects/preprocessing
rm -r Standalone_Wang2018_MFMem/stable_projects/registration
