#!/bin/sh
# Written by Lijun An and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd $CBIG_CODE_DIR/../
mkdir -p Standalone_An2024_DeepResBat/
rsync -axD $CBIG_CODE_DIR/stable_projects/harmonization/An2022_DeepResBat/ Standalone_An2024_DeepResBat/
