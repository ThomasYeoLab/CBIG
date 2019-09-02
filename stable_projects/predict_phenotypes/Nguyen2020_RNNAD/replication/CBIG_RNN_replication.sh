#!/bin/bash
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
ROOTDIR=`readlink -f $0 | xargs dirname`/../
export PYTHONPATH=$PYTHONPATH:$ROOTDIR

mkdir -p output

echo Generate 20 cross-validation folds
python -m cbig.Nguyen2020.gen_cv_fold \
    --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
    --features $ROOTDIR/data/features \
    --folds 20 \
    --outdir output

echo Create training and test data, model filling
python -m cbig.Nguyen2020.gen_cv_pickle \
    --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
    --features $ROOTDIR/data/features \
    --mask output/fold0_mask.csv \
    --strategy model \
    --batch_size 128 \
    --out output/test.f0.pkl

python -m cbig.Nguyen2020.gen_cv_pickle \
    --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
    --features $ROOTDIR/data/features \
    --mask output/fold1_mask.csv \
    --strategy model \
    --batch_size 128 \
    --out output/test.f1.pkl

echo Train MinimalRNN model, first fold
python -m cbig.Nguyen2020.train --verbose \
    --data output/test.f0.pkl \
    --i_drop 0.1 \
    --h_drop 0.4 \
    --h_size 128 \
    --nb_layers 2 \
    --epochs 100 --lr 0.001290666 --model MinRNN --weight_decay 1e-5 \
    --out output/model.f0.pt

echo Train MinimalRNN model, second fold
python -m cbig.Nguyen2020.train --verbose \
    --data output/test.f1.pkl \
    --i_drop 0.1 \
    --h_drop 0.4 \
    --h_size 128 \
    --nb_layers 2 \
    --epochs 100 --lr 0.001333218 --model MinRNN --weight_decay 1e-7 \
    --out output/model.f1.pt


echo Apply trained model on test set
python -m cbig.Nguyen2020.predict --checkpoint output/model.f0.pt --data output/test.f0.pkl \
    -o output/prediction_test.f0.csv
python -m cbig.Nguyen2020.predict --checkpoint output/model.f1.pt --data output/test.f1.pkl \
    -o output/prediction_test.f1.csv

echo; echo Evaluation prediction on test set, first fold
python -m cbig.Nguyen2020.evaluation --reference output/fold0_test.csv --prediction output/prediction_test.f0.csv
echo; echo Evaluation prediction on test set, second fold
python -m cbig.Nguyen2020.evaluation --reference output/fold1_test.csv --prediction output/prediction_test.f1.csv
