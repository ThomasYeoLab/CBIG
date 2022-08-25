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
import xgboost as xgb
from scipy.stats import sem
from config import global_config
from utils.misc import txt2list, list2txt, create_folder
from utils.normalization import ICV_norm
from evaluation.XGBoost.model import site_pred_model
from utils.metrics import site_prediction_metric


def train_XGBoost_args_parser():
    """
    Parameters for training XGBoost site prediction model
    """
    parser = argparse.ArgumentParser(prog='TrainXGBoostPredArgs')
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--output_path', type=str, default='/')
    parser.add_argument('--train_name', type=str, default='train')
    parser.add_argument('--val_name', type=str, default='val')
    parser.add_argument('--test_name', type=str, default='test')
    parser.add_argument('--save_suffix', type=str, default='/')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--norm', action='store_false', default=True)
    # model releated hyper-parameters
    parser.add_argument('--num_boost_rounds', type=int, default=100)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--threshold', type=float, default=0.01)
    # xgboost parameters
    params = {
        'booster': 'gbtree',
        'objective': 'binary:logistic',
        'eval_metric': 'error',
        'nthread': 1,
        'max_depth': 6,
        'subsample': 1.0
    }
    parser.add_argument('--params', type=dict, default=params)
    args, _ = parser.parse_known_args()
    return args


def grid_search(args, train_df, val_df, ROIs):
    """
    Perform grid search on validation set to get optimal hyper-parameters

    Args:
        args (tuple): Parameters
        train_df (class DataFrame): Dataframe for training
        val_df (class DataFrame): Dataframe for validation
        ROIs (list): Features
    """
    best_score, best_threshold, best_model = \
        site_pred_model(args, train_df, val_df, ROIs)
    # also outputs optimal hyper-parameters here is max_depth & subsample
    best_max_depth = 6
    best_subsample = 1.0
    for max_depth in (3, 4, 5, 6, 7):
        for subsample in (0.5, 0.6, 0.7, 0.8, 0.9, 1):
            args.params['max_depth'] = max_depth
            args.params['subsample'] = subsample
            score, threshold, model = \
                site_pred_model(args, train_df, val_df, ROIs)
            if score > best_score:
                best_score = score
                best_threshold = threshold
                best_model = model
                best_max_depth = max_depth
                best_subsample = subsample
    return best_model, best_threshold, best_max_depth, best_subsample


def train(args):
    """
    Train XGBoost model for 10 folds and save model & performance

    Args:
        args (tuple): Parameters
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
        if args.norm:
            # run ICV normalization
            norm_train_df = ICV_norm(train_df, ROIs)
            norm_val_df = ICV_norm(val_df, ROIs)
            norm_test_df = ICV_norm(test_df, ROIs)
        else:
            norm_train_df, norm_val_df, norm_test_df =\
                train_df, val_df, test_df
        # run grid search
        model, threshold, max_depth, subsamaple = \
            grid_search(args, norm_train_df, norm_val_df, ROIs)
        # save model and threshold
        model.save_model(
            os.path.join(fold_checkpoint_path,
                         'site_pred_model_' + args.save_suffix + '.json'))
        list2txt([threshold],
                 os.path.join(
                     fold_checkpoint_path,
                     'site_pred_threshold_' + args.save_suffix + '.txt'))
        list2txt([max_depth],
                 os.path.join(
                     fold_checkpoint_path,
                     'site_pred_max-depth_' + args.save_suffix + '.txt'))
        list2txt([subsamaple],
                 os.path.join(
                     fold_checkpoint_path,
                     'site_pred_subsample_' + args.save_suffix + '.txt'))
        # making prediction on testset and save prediction
        x_testset = norm_test_df[ROIs]
        y_testset = norm_test_df['SITE']
        xg_testset = xgb.DMatrix(x_testset, y_testset, feature_names=ROIs)
        pred = model.predict(xg_testset)
        acc_vec, acc = site_prediction_metric(norm_test_df['RID'].values, pred,
                                              y_testset.values, threshold)
        # save prediction
        np.save(
            os.path.join(fold_out_path,
                         'site_pred_test_' + args.save_suffix + '.npy'),
            acc_vec)
        logger.append(acc)
    acc_mean = np.mean(logger)
    acc_std = sem(logger)
    logger.append(str(acc_mean) + '_' + str(acc_std))
    # save logger
    list2txt(
        logger,
        os.path.join(args.output_path,
                     'site_pred_' + args.save_suffix + '.txt'))


if __name__ == '__main__':
    train(train_XGBoost_args_parser())
