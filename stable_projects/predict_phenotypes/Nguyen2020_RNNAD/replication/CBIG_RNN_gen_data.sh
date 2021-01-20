#!/bin/bash
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
ROOTDIR=`readlink -f "$0" | xargs dirname`/../
export PYTHONPATH=$PYTHONPATH:$ROOTDIR

mkdir -p data

echo Generate 20 cross-validation folds
python -m cbig.Nguyen2020.gen_cv_fold \
    --spreadsheet "$ROOTDIR"/data/TADPOLE_D1_D2.csv \
    --features "$ROOTDIR"/data/features \
    --folds 20 \
    --outdir data

echo Create training, validation and test data, model filling
for i in {0..19}; do
    python -m cbig.Nguyen2020.gen_cv_pickle \
        --spreadsheet "$ROOTDIR"/data/TADPOLE_D1_D2.csv \
        --features "$ROOTDIR"/data/features \
        --mask data/fold${i}_mask.csv \
        --strategy model \
        --batch_size 128 \
        --out data/val.f$i.pkl --validation

    python -m cbig.Nguyen2020.gen_cv_pickle \
        --spreadsheet "$ROOTDIR"/data/TADPOLE_D1_D2.csv \
        --features "$ROOTDIR"/data/features \
        --mask data/fold${i}_mask.csv \
        --strategy model \
        --batch_size 128 \
        --out data/test.f$i.pkl
done
