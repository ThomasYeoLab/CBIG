#!/usr/bin/env python
# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
import argparse
import os.path as path
import numpy as np
import pandas as pd

import cbig.Nguyen2020.misc as misc


def split_by_median_date(data, subjects):
    """
    Split timepoints in two halves, use first half to predict second half
    Args:
        data (Pandas data frame): input data
        subjects: list of subjects
    Return:
        first_half (ndarray): boolean mask, rows used as input
        second_half (ndarray): boolean mask, rows to predict
    """
    first_half = np.zeros(data.shape[0], int)
    second_half = np.zeros(data.shape[0], int)
    for rid in subjects:
        subj_mask = (data.RID == rid) & data.has_data
        median_date = np.sort(data.EXAMDATE[subj_mask])[subj_mask.sum() // 2]
        first_half[subj_mask & (data.EXAMDATE < median_date)] = 1
        second_half[subj_mask & (data.EXAMDATE >= median_date)] = 1
    return first_half, second_half


def gen_fold(data, nb_folds, outdir):
    """ Generate *nb_folds* cross-validation folds from *data """
    subjects = np.unique(data.RID)
    has_2tp = np.array([np.sum(data.RID == rid) >= 2 for rid in subjects])

    potential_targets = np.random.permutation(subjects[has_2tp])
    folds = np.array_split(potential_targets, nb_folds)

    leftover = [subjects[~has_2tp]]

    for test_fold in range(nb_folds):
        val_fold = (test_fold + 1) % nb_folds
        train_folds = [
            i for i in range(nb_folds) if (i != test_fold and i != val_fold)
        ]

        train_subj = np.concatenate(
            [folds[i] for i in train_folds] + leftover, axis=0)
        val_subj = folds[val_fold]
        test_subj = folds[test_fold]

        train_timepoints = (
            np.in1d(data.RID, train_subj) & data.has_data).astype(int)
        val_in_timepoints, val_out_timepoints = split_by_median_date(
            data, val_subj)
        test_in_timepoints, test_out_timepoints = split_by_median_date(
            data, test_subj)

        mask_frame = gen_mask_frame(data, train_timepoints, val_in_timepoints,
                                    test_in_timepoints)
        mask_frame.to_csv(
            path.join(outdir, 'fold%d_mask.csv' % test_fold), index=False)

        val_frame = gen_ref_frame(data, val_out_timepoints)
        val_frame.to_csv(
            path.join(outdir, 'fold%d_val.csv' % test_fold), index=False)

        test_frame = gen_ref_frame(data, test_out_timepoints)
        test_frame.to_csv(
            path.join(outdir, 'fold%d_test.csv' % test_fold), index=False)


def gen_mask_frame(data, train, val, test):
    """
    Create a frame with 3 masks:
        train: timepoints used for training model
        val: timepoints used for validation
        test: timepoints used for testing model
    """
    col = ['RID', 'EXAMDATE']
    ret = pd.DataFrame(data[col], index=range(train.shape[0]))
    ret['train'] = train
    ret['val'] = val
    ret['test'] = test

    return ret


def gen_ref_frame(data, test_timepoint_mask):
    """ Create reference frame which is used to evalute models' prediction """
    columns = [
        'RID', 'CognitiveAssessmentDate', 'Diagnosis', 'ADAS13', 'ScanDate'
    ]
    ret = pd.DataFrame(
        np.nan, index=range(len(test_timepoint_mask)), columns=columns)
    ret[columns] = data[['RID', 'EXAMDATE', 'DXCHANGE', 'ADAS13', 'EXAMDATE']]
    ret['Ventricles'] = data['Ventricles'] / data['ICV']
    ret = ret[test_timepoint_mask == 1]

    # map diagnosis from numeric categories back to labels
    mapping = {
        1: 'CN',
        7: 'CN',
        9: 'CN',
        2: 'MCI',
        4: 'MCI',
        8: 'MCI',
        3: 'AD',
        5: 'AD',
        6: 'AD'
    }
    ret.replace({'Diagnosis': mapping}, inplace=True)
    ret.reset_index(drop=True, inplace=True)

    return ret


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--spreadsheet', required=True)
    parser.add_argument('--features', required=True)
    parser.add_argument('--folds', type=int, required=True)
    parser.add_argument('--outdir', required=True)

    args = parser.parse_args()

    np.random.seed(args.seed)

    columns = ['RID', 'DXCHANGE', 'EXAMDATE']

    features = misc.load_feature(args.features)
    frame = pd.read_csv(
        args.spreadsheet,
        usecols=columns + features,
        converters=misc.CONVERTERS)
    frame['has_data'] = ~frame[features].isnull().apply(np.all, axis=1)
    gen_fold(frame, args.folds, args.outdir)


if __name__ == '__main__':
    main()
