#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd $CBIG_CODE_DIR/../
mkdir -p Standalone_An2022_gcVAE/
rsync -axD $CBIG_CODE_DIR/stable_projects/predict_phenotypes/An2022_gcVAE/ Standalone_An2022_gcVAE/
