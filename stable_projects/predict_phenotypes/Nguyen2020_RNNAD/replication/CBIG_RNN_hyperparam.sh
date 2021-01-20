#!/bin/bash
# Written by Minh Nguyen and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
ROOTDIR=`readlink -f $0 | xargs dirname`/../
export PYTHONPATH=$PYTHONPATH:$ROOTDIR

mkdir -p output

for i in {0..19}; do
    python -m cbig.Nguyen2020.hord_search \
        --executable 'python -m cbig.Nguyen2020.wrapper' \
        --data data/val.f$i.pkl \
        --epochs 100 --model MinRNN \
        --checkpoint output/model.f$i \
        --prediction output/pred_fold${i} \
        --reference data/fold${i}_val.csv \
        --maxeval 30 \
        --log_file output/summary_fold${i}.csv
done
