#!/bin/sh
#
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd $CBIG_CODE_DIR/../
mkdir -p Standalone_Nguyen2020_RNNAD/
rsync -axD $CBIG_CODE_DIR/stable_projects/predict_phenotypes/Nguyen2020_RNNAD/ Standalone_Nguyen2020_RNNAD/
