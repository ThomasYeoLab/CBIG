#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import numpy as np
from scipy import stats
from utils.misc import load_pkl, list2txt
from utils.metrics import subject_mae, subject_acc


def eval_goalDNN_args_parser():
    """
    Parameters for evaluating goalDNN's prediction

    Returns:
        tuple: Parameters
    """
    parser = argparse.ArgumentParser(prog='PredgoalDNNArgs')
    # parameters
    parser.add_argument('--save_path', type=str, default='/')
    parser.add_argument('--exp', type=str, default='sample_size')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--fold', type=int, default=0)
    parser.add_argument(
        '--harm_models', type=list, default=['ComBat', 'cVAE', 'gcVAE'])
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--testsets', type=list, default=[])

    args, _ = parser.parse_known_args()
    return args


def evaluation(pred_path, pred_names, nb_folds=10):
    """
    Evaluate the prediction performance of goalDNN

    Args:
        pred_path (str): Path for goalDNN prediction file
        pred_names (list): Models for performing harmonziation
        nb_folds (int, optional): Number of folds. Defaults to 10.
    """
    for pred_name in pred_names:
        mmse_mae_logger = []
        dx_acc_logger = []
        mmse_mae_array = []
        dx_acc_array = []
        for fold in range(nb_folds):
            fold_pred_path = os.path.join(pred_path, str(fold))
            pred_test = load_pkl(
                os.path.join(fold_pred_path,
                             'pred_test_' + pred_name + '.pkl'))
            # subject-level MMSE prediction MAE
            pred_test_mae_vector, pred_test_mae = \
                subject_mae(pred_test['MMSE']['Pred'],
                            pred_test['MMSE']['GT'],
                            pred_test['RID'].values,
                            'numpy')
            # subject-level Diagnosis prediction accuracy
            pred_test_acc_vector, pred_test_acc = \
                subject_acc(pred_test['DX']['Pred'],
                            pred_test['DX']['GT'],
                            pred_test['RID'].values,
                            'numpy')
            # logging results
            mmse_mae_logger.append(
                str(pred_test_mae) + '_' +
                str(stats.sem(pred_test_mae_vector)))
            dx_acc_logger.append(
                str(pred_test_acc) + '_' +
                str(stats.sem(pred_test_acc_vector)))
            mmse_mae_array.append(pred_test_mae_vector)
            dx_acc_array.append(pred_test_acc_vector)
        mmse_mae_array = np.array(mmse_mae_array)
        dx_acc_array = np.array(dx_acc_array)
        # firstly average across 10 folds, then average across subjects
        # mmse
        mmse_mean = np.mean(np.mean(mmse_mae_array, axis=0))
        mmse_sem = stats.sem(np.mean(mmse_mae_array, axis=0))
        mmse_mae_logger.append(str(mmse_mean) + '_' + str(mmse_sem))
        # dx
        dx_mean = np.mean(np.mean(dx_acc_array, axis=0))
        dx_sem = stats.sem(np.mean(dx_acc_array, axis=0))
        dx_acc_logger.append(str(dx_mean) + '_' + str(dx_sem))
        # save results
        list2txt(
            mmse_mae_logger,
            os.path.join(pred_path, 'MMSE_pred_result_' + pred_name + '.txt'))
        list2txt(
            dx_acc_logger,
            os.path.join(pred_path, 'DX_pred_result_' + pred_name + '.txt'))


def wrapper(args):
    """
    Wrapper function for evaluating goalDNN's prediction

    Args:
        args (tuple): Parameters
    """
    if args.exp == 'sample_size':
        pass
    else:
        args.harm_models = ['ComBat', 'ComBat4cov', 'cVAE', 'gcVAE']
    nonADNI_dataset = args.dataset_pair.split('-')[1]
    testsets = ['ADNI_unharm', nonADNI_dataset + '_unharm']
    for harm_model in args.harm_models:
        testsets.append(nonADNI_dataset + harm_model + '_harm')
    evaluation(args.save_path, testsets)


if __name__ == '__main__':
    wrapper(eval_goalDNN_args_parser())
