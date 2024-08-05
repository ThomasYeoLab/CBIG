#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import sem
from config import global_config
from model.XGB import XGBWrapper
from utils.misc import save_pkl, txt2list, list2txt, create_folder
from utils.normalization import ICV_norm
from utils.metrics import data_pred_metric


def dataset_pred_args_parser():
    """
    Parameters
    """
    parser = argparse.ArgumentParser(prog='DatasetPredArgs')
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--output_path', type=str, default='/')
    parser.add_argument('--train_name', type=str, default='train.csv')
    parser.add_argument('--val_name', type=str, default='val.csv')
    parser.add_argument('--test_name', type=str, default='test.csv')
    parser.add_argument('--sufix', type=str, default='/')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--num_boost_rounds', type=int, default=100)
    parser.add_argument('--seed', type=int, default=0)

    args, _ = parser.parse_known_args()
    return args


def dataset_pred_wrapper(args):
    """
    Wrapper function for predicting dataset
    """
    ROIs = txt2list(global_config.ROI_features_path)
    logger = []
    for fold in range(args.nb_folds):
        fold_input_path = os.path.join(args.data_path, str(fold))
        fold_checkpoint_path = os.path.join(args.checkpoint_path, str(fold))
        fold_out_path = os.path.join(args.output_path, str(fold))
        create_folder(fold_checkpoint_path)
        create_folder(fold_out_path)
        # read data
        train_df = pd.read_csv(
            os.path.join(fold_input_path, args.train_name + '.csv'))
        val_df = pd.read_csv(
            os.path.join(fold_input_path, args.val_name + '.csv'))
        test_df = pd.read_csv(
            os.path.join(fold_input_path, args.test_name + '.csv'))
        # run ICV normalization
        norm_train_df = ICV_norm(train_df, ROIs)
        norm_val_df = ICV_norm(val_df, ROIs)
        norm_test_df = ICV_norm(test_df, ROIs)

        # initialized XGBwrapper
        dataset_predictor = XGBWrapper(features=ROIs,
                                       target='SITE',
                                       task='site')
        # grid search
        best_model, best_hyper_params =\
            dataset_predictor.tune(norm_train_df, norm_val_df)
        # save
        checkpoint = dict()
        checkpoint['model'] = best_model
        checkpoint['best_hyper_params'] = best_hyper_params
        save_pkl(
            checkpoint,
            os.path.join(fold_checkpoint_path,
                         'dataset_pred_' + args.sufix + '.pkl'))
        # making prediction on testset
        test_pred = dataset_predictor.predict(best_model, norm_test_df)
        acc_vec, acc = data_pred_metric(norm_test_df['RID'].values, test_pred,
                                        norm_test_df['SITE'].values,
                                        best_hyper_params['threshold'])
        # save out acc_vec
        np.save(
            os.path.join(fold_out_path, 'dataset_pred_' + args.sufix + '.npy'),
            acc_vec)
        logger.append(acc)
    acc_mean = np.mean(logger)
    acc_std = sem(logger)
    logger.append(str(acc_mean) + '_' + str(acc_std))
    list2txt(
        logger,
        os.path.join(args.output_path, 'dataset_pred_' + args.sufix + '.txt'))


if __name__ == '__main__':
    dataset_pred_wrapper(dataset_pred_args_parser())
