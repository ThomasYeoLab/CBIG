#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
from __future__ import print_function
import argparse

import numpy as np

import cbig.Nguyen2020.misc as misc


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spreadsheet', required=True)
    parser.add_argument('--mask', '-d', required=True)
    parser.add_argument('--out', '-o', required=True)

    return parser.parse_args()


def last_value(observations, default):
    """ Get last observed values """
    has_data = ~observations.isnull()
    if sum(has_data) > 0:
        return observations[has_data].iloc[-1]

    return default


def to_categorical(last_status):
    """ Convert label to one-hot vector """
    if np.isnan(last_status):
        return np.ones(3) / 3

    ret = np.zeros(3)
    ret[int(last_status)] = 1
    return ret


def main(args):
    """ Constant baseline prediction """
    columns = ['RID', 'DX', 'ADAS13', 'Ventricles', 'ICV', 'EXAMDATE']
    frame = misc.load_table(args.spreadsheet, columns)
    # Ventricles volumes, normalised by intracranial volume
    frame['Vent'] = frame.Ventricles / frame.ICV

    train_mask, pred_mask, pred_mask_frame = misc.get_mask(
        args.mask, use_validation=False)
    assert len(frame) == len(train_mask)
    mean = frame[train_mask].mean(axis=0)

    _, start_dates = misc.get_baseline_prediction_start(pred_mask_frame)
    duration = 7 * 12

    subjects = np.unique(pred_mask_frame.RID)
    # 1. Clinical status forecasts
    dx = np.zeros([len(subjects), duration, 3])
    # 2. ADAS13 forecasts
    adas = np.zeros([len(subjects), duration, 3])
    # 3. Ventricles volume forecasts
    vent = np.zeros([len(subjects), duration, 3])

    for i, rid in enumerate(subjects):
        subj_frame = frame[(frame.RID == rid)
                           & pred_mask].sort_values('EXAMDATE')

        dx[i] = to_categorical(last_value(subj_frame.DX, np.nan))
        adas[i] = misc.add_ci_col(
            last_value(subj_frame.ADAS13, mean.ADAS13), 1, 0, 85)
        vent[i] = misc.add_ci_col(
            last_value(subj_frame.Vent, mean.Vent), 5e-4, 0., 1.)

    prediction = {'subjects': subjects}
    prediction['dates'] = misc.make_date_col(
        [start_dates[rid] for rid in subjects], duration)
    prediction['DX'] = dx
    prediction['ADAS13'] = adas
    prediction['Ventricles'] = vent

    misc.build_pred_frame(prediction, args.out)


if __name__ == '__main__':
    main(get_args())
