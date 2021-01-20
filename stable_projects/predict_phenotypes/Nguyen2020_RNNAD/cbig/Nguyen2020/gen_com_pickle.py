#!/usr/bin/env python
# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
from __future__ import print_function, division
import argparse
import pickle

import cbig.Nguyen2020.dataloader as dataloader
import cbig.Nguyen2020.misc as misc


def get_data(args, fields, mean=None, stds=None):
    """
    Generate training/test data batches and save as pickle file
    *args* specify
        mask: path to mask file
        strategy: filling strategy
        spreadsheet: path to TADPOLE data
        start_date: date when prediction starts
        nb_years: number of years to predict
        batch_size: batch size
        out: path to save pickle file
    """

    ret = {}
    train_mask, pred_mask, pred_mask_frame = misc.get_mask(args.mask, False)
    ret['baseline'], _ = misc.get_baseline_prediction_start(pred_mask_frame)
    prediction_start = misc.str2date(args.start_date)
    ret['pred_start'] = {rid: prediction_start for rid in ret['baseline']}
    ret['duration'] = args.nb_years * 12

    columns = ['RID', 'Month_bl', 'DX'] + fields
    frame = misc.load_table(args.spreadsheet, columns)

    tf = frame.loc[train_mask, fields]
    ret['mean'] = tf.mean() if mean is None else mean
    ret['stds'] = tf.std() if stds is None else stds
    ret['VentICVstd'] = (tf['Ventricles'] / tf['ICV']).std()

    frame[fields] = (frame[fields] - ret['mean']) / ret['stds']

    default_val = {f: 0. for f in fields}
    default_val['DX'] = 0.

    data = dataloader.extract(frame[train_mask], args.strategy, fields,
                              default_val)
    ret['train'] = dataloader.Random(data, args.batch_size, fields)

    data = dataloader.extract(frame[pred_mask], args.strategy, fields,
                              default_val)
    ret['test'] = dataloader.Sorted(data, 1, fields)

    print('train', len(ret['train'].subjects), 'subjects')
    print('test', len(ret['test'].subjects), 'subjects')
    print(len(fields), 'features')

    return ret


def main(args):
    if args.feat_stat is not None:
        with open(args.feat_stat, 'rb') as fhandler:
            mean, stds = pickle.load(fhandler)
        fields = list(mean.index)
    else:
        mean, stds = None, None
        fields = misc.load_feature(args.features)
    with open(args.out, 'wb') as fhandler:
        pickle.dump(get_data(args, fields, mean, stds), fhandler)


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mask', '-m', required=True)
    parser.add_argument('--strategy', '-s', required=True)
    parser.add_argument('--spreadsheet', '-t', required=True)
    parser.add_argument('--features')
    parser.add_argument('--feat_stat')
    parser.add_argument('--start_date', default='2018-01-01')
    parser.add_argument('--nb_years', type=int, default=5)
    parser.add_argument('--batch_size', type=int, required=True)
    parser.add_argument('--out', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    main(get_args())
