#!/bin/bash
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
ROOTDIR=`readlink -f $0 | xargs dirname`/../
export PYTHONPATH=$PYTHONPATH:$ROOTDIR

mkdir -p output

echo Generate 10 cross-validation folds
python -m cbig.Nguyen2020.gen_cv_fold \
    --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
    --features $ROOTDIR/data/features \
    --folds 10 \
    --outdir output

echo Create training and validation data using the first fold, model filling
python -m cbig.Nguyen2020.gen_cv_pickle \
    --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
    --features $ROOTDIR/data/features \
    --mask output/fold0_mask.csv \
    --strategy model \
    --batch_size 128 \
    --out output/val.pkl --validation

echo Train MinimalRNN model using the training set of the first fold
python -m cbig.Nguyen2020.train --verbose \
    --data output/val.pkl \
    --i_drop 0.1 \
    --h_drop 0.1 \
    --h_size 512 \
    --epochs 100 --lr 5e-4 --model MinRNN --weight_decay 5e-7 \
    --out output/model.pt

echo Apply trained model on validation set
python -m cbig.Nguyen2020.predict --checkpoint output/model.pt --data output/val.pkl -o output/prediction_val.csv

echo Evaluation prediction on validation set
python -m cbig.Nguyen2020.evaluation --reference output/fold0_val.csv --prediction output/prediction_val.csv

echo Create training and test data using the first fold, model filling
python -m cbig.Nguyen2020.gen_cv_pickle \
    --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
    --features $ROOTDIR/data/features \
    --mask output/fold0_mask.csv \
    --strategy model \
    --batch_size 128 \
    --out output/test.pkl

echo Apply trained model on test set
python -m cbig.Nguyen2020.predict --checkpoint output/model.pt --data output/test.pkl -o output/prediction_test.csv

echo Evaluation prediction on test set
python -m cbig.Nguyen2020.evaluation --reference output/fold0_test.csv --prediction output/prediction_test.csv

echo Constant baseline
python -m cbig.Nguyen2020.baseline_constant \
     --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
     --mask output/fold0_mask.csv \
     --out output/constant_prediction.csv

python -m cbig.Nguyen2020.evaluation --reference output/fold0_test.csv --prediction output/constant_prediction.csv

echo SVM baseline
python -m cbig.Nguyen2020.baseline_svm \
     --spreadsheet $ROOTDIR/data/TADPOLE_D1_D2.csv \
     --mask output/fold0_mask.csv \
     --features $ROOTDIR/data/features \
     --agamma .01 --vgamma .01 \
     --out output/svm_prediction.csv

python -m cbig.Nguyen2020.evaluation --reference output/fold0_test.csv --prediction output/svm_prediction.csv
