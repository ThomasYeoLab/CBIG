#!/bin/bash
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
ROOTDIR=`readlink -f $0 | xargs dirname`/../
export PYTHONPATH=$PYTHONPATH:$ROOTDIR

echo Create mask for D2 subjects
python create_mask.py ../data/TADPOLE_D1_D2.csv tadpole_mask.csv

echo Generate test data
python -m cbig.Nguyen2020.gen_com_pickle \
    --mask tadpole_mask.csv \
    --strategy model \
    --spreadsheet ../data/TADPOLE_D1_D2.csv \
    --feat_stat save.seed0/feat_stats.pkl \
    --batch_size 128 \
    --out tadpole_test.pkl

echo Predict using 1st set of pre-trained weights
python -m cbig.Nguyen2020.predict \
    --checkpoint save.seed0 \
    --data tadpole_test.pkl \
    -o prediction0.csv

echo Predict using 2nd set of pre-trained weights
python -m cbig.Nguyen2020.predict \
    --checkpoint save.seed1 \
    --data tadpole_test.pkl \
    -o prediction1.csv

echo Predict using 3rd set of pre-trained weights
python -m cbig.Nguyen2020.predict \
    --checkpoint save.seed2 \
    --data tadpole_test.pkl \
    -o prediction2.csv

echo Predict using 4th set of pre-trained weights
python -m cbig.Nguyen2020.predict \
    --checkpoint save.seed3 \
    --data tadpole_test.pkl \
    -o prediction3.csv

python avg_prediction.py -f prediction{0,1,2,3}.csv -s D4Live
