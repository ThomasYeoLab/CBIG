#!/usr/bin/env python
# Written by Minh Nguyen and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
from __future__ import print_function, division
import argparse
import itertools

import numpy as np
import cbig.Nguyen2020.misc as misc


def a_value(probabilities, zero_label=0, one_label=1):
    """
    Approximates the AUC by the method described in Hand and Till 2001,
    equation 3.
    NB: The class labels should be in the set [0,n-1] where n = # of classes.
    The class probability should be at the index of its label in the
    probability list.
    I.e. With 3 classes the labels should be 0, 1, 2. The class probability
    for class '1' will be found in index 1 in the class probability list
    wrapped inside the zipped list with the labels.
    Args:
        probabilities (list): A zipped list of the labels and the
            class probabilities in the form (m = # data instances):
             [(label1, [p(x1c1), p(x1c2), ... p(x1cn)]),
              (label2, [p(x2c1), p(x2c2), ... p(x2cn)])
                             ...
              (labelm, [p(xmc1), p(xmc2), ... (pxmcn)])
             ]
        zero_label (optional, int): The label to use as the class '0'.
            Must be an integer, see above for details.
        one_label (optional, int): The label to use as the class '1'.
            Must be an integer, see above for details.
    Returns:
        The A-value as a floating point.
    Source of script: https://gist.github.com/stulacy/672114792371dc13b247
    """
    # Obtain a list of the probabilities for the specified zero label class
    expanded_points = []
    for instance in probabilities:
        if instance[0] == zero_label or instance[0] == one_label:
            expanded_points.append((instance[0], instance[1][zero_label]))
    sorted_ranks = sorted(expanded_points, key=lambda x: x[1])

    n0, n1, sum_ranks = 0, 0, 0
    # Iterate through ranks and increment counters for overall count and ranks
    for index, point in enumerate(sorted_ranks):
        if point[0] == zero_label:
            n0 += 1
            sum_ranks += index + 1  # Add 1 as ranks are one-based
        elif point[0] == one_label:
            n1 += 1
        else:
            pass  # Not interested in this class

    return (sum_ranks - (n0 * (n0 + 1) / 2.0)) / float(n0 * n1)  # Eqn 3


def MAUC(data, no_classes):
    """
    Calculates the MAUC over a set of multi-class probabilities and
    their labels. This is equation 7 in Hand and Till's 2001 paper.
    NB: The class labels should be in the set [0,n-1] where n = # of classes.
    The class probability should be at the index of its label in the
    probability list.
    I.e. With 3 classes the labels should be 0, 1, 2. The class probability
    for class '1' will be found in index 1 in the class probability list
    wrapped inside the zipped list with the labels.
    Args:
        data (list): A zipped list (NOT A GENERATOR) of the labels and the
            class probabilities in the form (m = # data instances):
             [(label1, [p(x1c1), p(x1c2), ... p(x1cn)]),
              (label2, [p(x2c1), p(x2c2), ... p(x2cn)])
                             ...
              (labelm, [p(xmc1), p(xmc2), ... (pxmcn)])
             ]
        no_classes (int): The number of classes in the dataset.
    Returns:
        The MAUC as a floating point value.
    Source of script: https://gist.github.com/stulacy/672114792371dc13b247
    """
    # Find all pairwise comparisons of labels
    class_pairs = [x for x in itertools.combinations(range(no_classes), 2)]

    # Have to take average of A value with both classes acting as label 0 as
    # this gives different outputs for more than 2 classes
    sum_avals = 0
    for pairing in class_pairs:
        sum_avals += (
            a_value(data, zero_label=pairing[0], one_label=pairing[1]) +
            a_value(data, zero_label=pairing[1], one_label=pairing[0])) / 2.0

    return sum_avals * (2 / float(no_classes * (no_classes - 1)))  # Eqn 7


def calcBCA(estimLabels, trueLabels, no_classes):
    """
    Calculates the balanced class accuracy (BCA)
    Args:
        estimLabels (ndarray): predicted classes
        trueLabels (ndarray): ground truth classes
        no_classes (int): The number of classes in the dataset.
    Returns:
        BCA value
    """
    bcaAll = []
    for c0 in range(no_classes):
        # c0 can be either CTL, MCI or AD

        # one example when c0=CTL
        # TP - label was estimated as CTL, and the true label was also CTL
        # FP - label was estimated as CTL, but the true label was not CTL
        TP = np.sum((estimLabels == c0) & (trueLabels == c0))
        TN = np.sum((estimLabels != c0) & (trueLabels != c0))
        FP = np.sum((estimLabels == c0) & (trueLabels != c0))
        FN = np.sum((estimLabels != c0) & (trueLabels == c0))

        # sometimes the sensitivity of specificity can be NaN, if the user
        # doesn't forecast one of the classes.
        # In this case we assume a default value for sensitivity/specificity
        if (TP + FN) == 0:
            sensitivity = 0.5
        else:
            sensitivity = (1. * TP) / (TP + FN)

        if (TN + FP) == 0:
            specificity = 0.5
        else:
            specificity = (1. * TN) / (TN + FP)

        bcaCurr = 0.5 * (sensitivity + specificity)
        bcaAll += [bcaCurr]

    return np.mean(bcaAll)


def nearest(series, target):
    """ Return index in *series* with value closest to *target* """
    return (series - target).abs().idxmin()


def mask(pred, true):
    """ Drop entries without ground truth data (i.e. NaN in *true*) """
    try:
        index = ~np.isnan(true)
    except Exception:
        print('true', true)
        print('pred', pred)
        raise
    ret = pred[index], true[index]
    assert ret[0].shape[0] == ret[0].shape[0]
    return ret


def parse_data(_ref_frame, _pred_frame):
    """ Match ground truth timepoints to closest predicted timepoints """
    true_label_and_prob = []
    pred_diag = np.full(len(_ref_frame), -1, dtype=int)
    pred_adas = np.full(len(_ref_frame), -1, dtype=float)
    pred_vent = np.full(len(_ref_frame), -1, dtype=float)

    for i in range(len(_ref_frame)):
        cur_row = _ref_frame.iloc[i]
        subj_data = _pred_frame[_pred_frame.RID == cur_row.RID].reset_index(
            drop=True)
        dates = subj_data['Forecast Date']

        matched_row = subj_data.iloc[nearest(dates,
                                             cur_row.CognitiveAssessmentDate)]
        prob = matched_row[[
            'CN relative probability', 'MCI relative probability',
            'AD relative probability'
        ]].values
        pred_diag[i] = np.argmax(prob)
        pred_adas[i] = matched_row['ADAS13']

        # for the mri scan find the forecast closest to the scan date,
        # which might be different from the cognitive assessment date
        pred_vent[i] = subj_data.iloc[nearest(
            dates, cur_row.ScanDate)]['Ventricles_ICV']

        if not np.isnan(cur_row.Diagnosis):
            true_label_and_prob += [(cur_row.Diagnosis, prob)]

    pred_diag, true_diag = mask(pred_diag, _ref_frame.Diagnosis)
    pred_adas, true_adas = mask(pred_adas, _ref_frame.ADAS13)
    pred_vent, true_vent = mask(pred_vent, _ref_frame.Ventricles)

    return true_label_and_prob, pred_diag, pred_adas, pred_vent, \
        true_diag, true_adas, true_vent


def is_date_column(col):
    """ Is the column of type datetime """
    return np.issubdtype(col.dtype, np.datetime64)


def eval_submission(ref_frame, pred_frame):
    """ Evaluate mAUC, BCA, ADAS13 MAE, and ventricles MAE """
    assert is_date_column(ref_frame['CognitiveAssessmentDate'])
    assert is_date_column(ref_frame['ScanDate'])
    assert is_date_column(pred_frame['Forecast Date'])

    true_labels_and_prob, p_diag, p_adas, p_vent, t_diag, t_adas, t_vent = \
        parse_data(ref_frame, pred_frame)

    try:
        mauc = MAUC(true_labels_and_prob, no_classes=3)
    except ZeroDivisionError:
        mauc = float('NaN')

    bca = calcBCA(p_diag, t_diag.astype(int), no_classes=3)

    adas = np.mean(np.abs(p_adas - t_adas))
    vent = np.mean(np.abs(p_vent - t_vent))

    return {'mAUC': mauc, 'bca': bca, 'adasMAE': adas, 'ventsMAE': vent}


def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', '-r', required=True)
    parser.add_argument('--prediction', '-p', required=True)
    return parser


def main(args):
    result = eval_submission(
        misc.read_csv(args.reference), misc.read_csv(args.prediction))

    print('########### Metrics for clinical status ##################')
    print('mAUC', result['mAUC'], 'bca', result['bca'])
    print('########### Mean Absolute Error (MAE) ##################')
    print('adasMAE', result['adasMAE'], 'ventsMAE', result['ventsMAE'])


if __name__ == "__main__":
    main(get_arg_parser().parse_args())
