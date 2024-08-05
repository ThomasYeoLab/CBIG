#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import copy
import argparse
import pandas as pd
from config import global_config
from utils.imputation import forward_filling
from utils.normalization import z_norm
from utils.misc import create_folder, txt2list, save_pkl


def deepharm_input_gen_args_parser():
    """
    Parameters for generating input for VAE models
    """
    parser = argparse.ArgumentParser(prog='VAEInputArgs')
    parser.add_argument('--input_path', type=str, default='/')
    parser.add_argument('--harm_input_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--train_name',
                        type=str,
                        default='unmatch2match_train')
    parser.add_argument('--val_name', type=str, default='unmatch2match_val')
    parser.add_argument('--test_name', type=str, default='unmatch2match_test')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--demean', action='store_false', default=True)

    args, _ = parser.parse_known_args()
    return args


def deepharm_input_gen(args):
    """
    Generate input for VAE models
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    numerical_features = ROIs + ['AGE', 'SITE', 'SEX', 'MMSE', 'DX']
    categorical_features = ['SEX', 'DX', 'SITE']
    for fold in range(args.nb_folds):
        fold_input_path = os.path.join(args.input_path, args.dataset_pair,
                                       str(fold))
        fold_output_path = os.path.join(args.harm_input_path,
                                        args.dataset_pair, str(fold))
        create_folder(fold_output_path)
        # Train
        train_file_name = args.train_name + '.csv'
        train_df = pd.read_csv(os.path.join(fold_input_path, train_file_name))
        train_df = train_df[cols]
        train_data = copy.deepcopy(train_df)
        train = dict()
        train['mean'] = train_data[numerical_features].mean()
        train['std'] = train_data[numerical_features].std()
        train_data = forward_filling(train_data, 'MMSE')
        train_data = forward_filling(train_data, 'DX')
        train_data.fillna(value=0, inplace=True)
        train['mean'][categorical_features] = 0.
        train['std'][categorical_features] = 1.
        if not args.demean:
            train['mean'][ROIs] = 0
        train['data'] = z_norm(train_data[numerical_features], train['mean'],
                               train['std'])
        train['RID'] = train_data[['RID']]
        save_pkl(train, os.path.join(fold_output_path,
                                     args.train_name + '.pkl'))
        # Validation
        val_file_name = args.val_name + '.csv'
        val_df = pd.read_csv(os.path.join(fold_input_path, val_file_name))
        val = dict()
        val['mean'] = train['mean']
        val['std'] = train['std']
        val_data = copy.deepcopy(val_df)  # unmatched for tunning
        val_data = forward_filling(val_data, 'MMSE')
        val_data = forward_filling(val_data, 'DX')
        val_data.fillna(value=0, inplace=True)
        val['data'] = z_norm(val_data[numerical_features], val['mean'],
                             val['std'])
        val['RID'] = val_data[['RID']]
        save_pkl(val, os.path.join(fold_output_path, args.val_name + '.pkl'))
        # Test
        test_file_name = args.test_name + '.csv'
        test_df = pd.read_csv(os.path.join(fold_input_path, test_file_name))
        test = dict()
        test['mean'] = train['mean']
        test['std'] = train['std']
        test_data = copy.deepcopy(test_df)  # unmatched for tunning
        test_data = forward_filling(test_data, 'MMSE')
        test_data = forward_filling(test_data, 'DX')
        test_data.fillna(value=0, inplace=True)
        test['data'] = z_norm(test_data[numerical_features], test['mean'],
                              test['std'])
        test['RID'] = test_data[['RID']]
        save_pkl(test, os.path.join(fold_output_path, args.test_name + '.pkl'))


if __name__ == '__main__':
    deepharm_input_gen(deepharm_input_gen_args_parser())
