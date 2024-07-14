#!/bin/sh
# Written by NAREN WULAN and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd $CBIG_CODE_DIR/../
mkdir -p Standalone_Naren2024_MMT1/
rsync -axD $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Naren2024_MMT1/ Standalone_Naren2024_MMT1/
