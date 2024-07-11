#!/bin/sh
#
# Written by Pansheng Chen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ~/storage
rsync -a --exclude .git $CBIG_CODE_DIR/* Standalone_Chen2024_MMM

# remove useless stable projects
rm -r Standalone_Chen2024_MMM/stable_projects/brain_parcellation
rm -r Standalone_CHen2024_MMM/stable_projects/disorder_subtypes
rm -r Standalone_Chen2024_MMM/stable_projects/fMRI_dynamics
rm -r Standalone_Chen2024_MMM/stable_projects/meta-analysis
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/An2022_gc_VAE
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/ChenOoi2023_ICCW
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/ChenTam2022_TRBPC
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/He2019_KRDNN
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/He2022_MM
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/Kong2023_GradPar
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/Nguyen2020_RNNAD
rm -r Standalone_Chen2024_MMM/stable_projects/predict_phenotypes/Ooi2022_MMP
rm -r Standalone_Chen2024_MMM/stable_projects/preprocessing/CBIG_fMRI_Preproc2016
rm -r Standalone_Chen2024_MMM/stable_projects/preprocessing/CBIG2022_DiffProc
rm -r Standalone_Chen2024_MMM/stable_projects/preprocessing/Orban2020_tod
rm -r Standalone_Chen2024_MMM/stable_projects/registration
