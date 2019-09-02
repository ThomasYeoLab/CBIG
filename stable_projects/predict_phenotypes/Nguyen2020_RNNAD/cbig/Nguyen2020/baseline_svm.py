#!/usr/bin/env python
# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
from __future__ import print_function

import argparse
from dateutil.relativedelta import relativedelta

from joblib import Parallel, delayed
import pandas as pd
import numpy as np
import sklearn.svm as svm

import cbig.Nguyen2020.misc as misc

WORKERS = Parallel(n_jobs=8)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spreadsheet', required=True)
    parser.add_argument('--mask', required=True)
    parser.add_argument('--features', required=True)

    parser.add_argument('--timepoints', '-t', type=int, default=1)

    parser.add_argument('--dkernel', choices=['linear', 'rbf'], default='rbf')
    parser.add_argument('--dc', type=float, default=1e0)
    parser.add_argument('--dgamma', type=float, default=1e-1)

    parser.add_argument('--akernel', choices=['linear', 'rbf'], default='rbf')
    parser.add_argument('--ac', type=float, default=1e0)
    parser.add_argument('--agamma', type=float, default=1e0)
    parser.add_argument('--aeps', type=float, default=1e-1)

    parser.add_argument('--vkernel', choices=['linear', 'rbf'], default='rbf')
    parser.add_argument('--vc', type=float, default=1e0)
    parser.add_argument('--vgamma', type=float, default=1e0)
    parser.add_argument('--veps', type=float, default=1e-1)

    parser.add_argument('--out', '-o', required=True)

    return parser.parse_args()


def get_offset(start, interval, count):
    """ Return *count* offsets which are *interval* apart """
    return start + np.arange(count) * interval


def interp(frame, idx, **args):
    """ Get interpolated values at indices """
    unique_idx = [i for i in idx if i not in frame.index]
    aug_frame = frame.append(pd.DataFrame(index=unique_idx), sort=True)
    in_frame = aug_frame.sort_index().interpolate('index', **args)

    return in_frame.loc[idx, frame.columns]


def get_subj_traindata(frame, offsets_from_tgt, fields):
    """ Get training data for each subjects """
    frames = [
        interp(frame, frame.index - idx, limit_direction='backward')
        for idx in offsets_from_tgt
    ]
    in_blk = [x.values for x in frames]

    mask = np.all(
        np.vstack([(x.index >= 0) & misc.has_data_mask(x) for x in frames]),
        axis=0)

    return np.hstack(in_blk)[mask], {
        f: frame.loc[mask, f].values
        for f in fields
    }


def get_traindata(threadpool, data, offsets_from_tgt, fields):
    """ Get training data """
    jobs = [
        delayed(get_subj_traindata)(sf, offsets_from_tgt, fields)
        for sf in data.values()
    ]
    input_, output_ = zip(*threadpool(jobs))

    input_ = np.concatenate(input_)
    output_ = {f: np.concatenate([x[f] for x in output_]) for f in fields}

    return input_, output_


def fit(tags, input_, output, **params):
    """ Train individual model """
    model_class = params['class']
    model_params = dict(params)
    model_params.pop('class')
    model = model_class(**model_params)

    valid = ~np.isnan(output)
    if np.any(valid):
        model.fit(input_[valid], output[valid])
    return model, tags


def get_test_input(data, subjects, offsets_from_inp):
    """ Get prediction input """
    input_ = []
    for subj in subjects:
        sf = data[subj]
        frame = interp(
            sf, max(sf.index) - offsets_from_inp, limit_area='inside')
        input_.append(np.hstack([x for x in frame.values]))

    return np.vstack(input_)


def predict_helper(model, data):
    """ Return SVM/SVR prediction """
    if isinstance(model, svm.SVC):
        return model.predict_proba(data)
    elif isinstance(model, svm.SVR):
        return model.predict(data)
    raise TypeError('model')


def keypoints(frame, interval):
    """ Return indices starting from last indices that are *interval* apart """
    return max(frame.index) - min(frame.index) // interval


def select(timepoints, prediction):
    """ Select the prediction that uses *timepoints* number of timepoints """
    tp_indices = timepoints - 1
    subj_indices = np.arange(len(tp_indices))
    return prediction[tp_indices, subj_indices]


def group_predict(data, max_tp, interval, model_groups):
    """ Predict using multiple models at multiple timepoints """
    assert isinstance(model_groups, dict)
    assert np.all([max_tp == len(models) for models in model_groups.values()])

    subjects = sorted(data.keys())

    predictions = {k: [] for k in model_groups}
    for i in range(max_tp):
        input_ = get_test_input(data, subjects, get_offset(0, interval, i + 1))
        input_[np.isnan(input_)] = 0.

        for j, models in model_groups.items():
            predictions[j].append(predict_helper(models[i], input_))

    timepoints = [keypoints(data[s], interval) for s in subjects]
    timepoints = np.array(timepoints, dtype=int).clip(max=max_tp)

    return {k: select(timepoints, np.array(p)) for k, p in predictions.items()}


def add_ci_col(values, c_interval, lowerbound, upperbound):
    """" Add connfidence intervals """
    val = np.expand_dims(values, -1)
    lo_ci = val - c_interval
    hi_ci = val + c_interval
    out = np.concatenate([val, lo_ci, hi_ci], axis=-1)
    return np.clip(out, lowerbound, upperbound)


def fit_models(data, max_tp, interval, nb_predictions, model_params):
    """ Train SVM/SVR models """
    targets = model_params.keys()
    models = {k: [] for k in targets}
    for i in range(nb_predictions):
        for j in model_params:
            models[j].append([None] * max_tp)

    jobs = []
    for i in range(nb_predictions):
        for j in range(max_tp):
            indices = get_offset(interval * (i + 1), interval, j + 1)
            in_, out_ = get_traindata(WORKERS, data, indices, targets)
            in_[np.isnan(in_)] = 0.

            for k, params in model_params.items():
                tag = (k, i, j)
                jobs.append(delayed(fit)(tag, in_, out_[k], **params))

    for model, (k, i, j) in WORKERS(jobs):
        models[k][i][j] = model

    return models


def extract(out_matrix, boundaries):
    """ Extract matrices within boundaries """
    assert out_matrix.shape[0] == len(boundaries)
    ret = []
    for i, (start, end) in enumerate(boundaries):
        ret.append(out_matrix[i:i + 1, start:end])
    return np.vstack(ret)


def get_boundaries(data, pred_frame):
    """
    Get start and end indices of prediction and
    dates corresponding to this range
    """
    baselines, starts = misc.get_baseline_prediction_start(pred_frame)
    duration = 7 * 12

    boundaries, dates = [], []
    interval_range = [relativedelta(months=i) for i in range(duration)]
    for subj in sorted(data.keys()):
        start = starts[subj]
        idx_diff = misc.month_between(start, baselines[subj])
        start_idx = idx_diff - max(data[subj].index)
        assert start_idx >= 0
        boundaries.append((start_idx, start_idx + duration))

        dates.append([start + d for d in interval_range])

    return boundaries, dates


def interp_output(diag, adas, vent, indices, boundaries):
    """ Linear interpolate output. Extrapolate after last prediction """
    idx = np.arange(150)

    def get_segment(in_frame):
        unique_idx = [i for i in idx if i not in in_frame.index]
        aug_frame = in_frame.append(pd.DataFrame(index=unique_idx), sort=True)
        out_frame = aug_frame.sort_index().interpolate(
            'index', limit_direction='both')

        return extract(out_frame.loc[idx, in_frame.columns].values.T,
                       boundaries)

    out_diag = []
    diag_cat = np.array(diag)
    for i in range(diag_cat.shape[-1]):
        out_diag.append(
            get_segment(pd.DataFrame(diag_cat[:, :, i], index=indices)))

    out_adas = get_segment(pd.DataFrame(np.vstack(adas), index=indices))
    out_vent = get_segment(pd.DataFrame(np.vstack(vent), index=indices))

    ret = {'DX': np.concatenate([x[:, :, None] for x in out_diag], axis=-1)}
    ret['ADAS13'] = add_ci_col(out_adas, 1., 0, 85)
    ret['Ventricles'] = add_ci_col(out_vent, 0.05, 0, 1)

    return ret


def pack(dict_of_list):
    """ Return dictionary of models with key is a tuple of (str, int) """
    ret = {}
    for k in dict_of_list:
        for i, x in enumerate(dict_of_list[k]):
            ret[(k, i)] = x
    return ret


def unpack(dict_of_list, dict_):
    """
    Unpack predictions into a dictionary.
    Dictionary keys are keys of *dict_*
    """
    ret = {k: [None] * len(v) for k, v in dict_of_list.items()}
    for k, p in dict_.items():
        ret[k[0]][k[1]] = p
    return ret


def get_model_params(args):
    """ Extract model parameters from command line arguments """
    dx_params = {'class': svm.SVC, 'max_iter': 100000, 'probability': True}
    dx_params['kernel'] = args.dkernel
    dx_params['C'] = args.dc
    dx_params['gamma'] = args.dgamma

    vent_params = {'class': svm.SVR, 'max_iter': 100000}
    vent_params['kernel'] = args.vkernel
    vent_params['C'] = args.vc
    vent_params['gamma'] = args.vgamma
    vent_params['epsilon'] = args.veps

    adas_params = {'class': svm.SVR, 'max_iter': 100000}
    adas_params['kernel'] = args.akernel
    adas_params['C'] = args.ac
    adas_params['gamma'] = args.agamma
    adas_params['epsilon'] = args.aeps

    return {'DX': dx_params, 'ADAS13': adas_params, 'Vent': vent_params}


def get_last_value(data, subjects):
    """ Get last observed values """
    last_values = []
    for id_ in subjects:
        last_tp = data[id_].interpolate(
            'index', limit_direction='forward').iloc[-1]
        last_values.append(last_tp)
    return pd.concat(last_values)


def main(args):
    """ SVM baseline """
    np.random.seed(0)
    features = misc.load_feature(args.features)

    frame = misc.load_table(args.spreadsheet,
                            ['RID', 'DX', 'Month_bl'] + features)
    frame['Vent'] = frame.Ventricles / frame.ICV
    features.append('Vent')

    train_tp, pred_tp, pred_frame = misc.get_mask(
        args.mask, use_validation=False)

    mean = frame.loc[train_tp, features].mean()
    std = frame.loc[train_tp, features].std()

    train_df = frame[train_tp].copy()
    test_df = frame[pred_tp].copy()

    train_df[features] = (train_df[features] - mean) / std
    test_df[features] = (test_df[features] - mean) / std

    # Training
    interval, nb_keypoints = 6, 10
    model_groups = fit_models(
        misc.get_data_dict(train_df, features), args.timepoints, interval,
        nb_keypoints, get_model_params(args))

    # Prediction
    test_data = misc.get_data_dict(test_df, features)
    pred = {'subjects': np.array(sorted(test_data.keys()))}
    last_values = get_last_value(test_data, pred['subjects'])

    predictions = group_predict(test_data, args.timepoints, interval,
                                pack(model_groups))
    predictions = unpack(model_groups, predictions)
    dx_keypoints = [misc.to_categorical(last_values.DX, 3)] + predictions['DX']
    adas_keypoints = [last_values.ADAS13] + predictions['ADAS13']
    vent_keypoints = [last_values.Vent] + predictions['Vent']

    adas_keypoints = [x * std.ADAS13 + mean.ADAS13 for x in adas_keypoints]
    vent_keypoints = [x * std.Vent + mean.Vent for x in vent_keypoints]

    boundaries, pred['dates'] = get_boundaries(test_data, pred_frame)

    pred.update(
        interp_output(dx_keypoints, adas_keypoints, vent_keypoints,
                      get_offset(0, interval, len(adas_keypoints)),
                      boundaries))

    misc.build_pred_frame(pred, args.out)


if __name__ == '__main__':
    main(get_args())
