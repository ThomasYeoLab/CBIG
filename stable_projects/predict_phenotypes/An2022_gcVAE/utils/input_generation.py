#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import copy
import argparse
import warnings
import pandas as pd
from config import global_config
from utils.normalization import z_norm
from utils.imputation import forward_filling
from utils.misc\
    import create_folder, txt2list, save_pkl, load_pkl, csvs2df

warnings.filterwarnings('ignore')


def input_gen_args_parser():
    """
    Parameters for generating input

    Returns:
        _type_: _description_
    """
    parser = argparse.ArgumentParser(prog='GenInputArgs')
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--exp', type=str, default='unmatch2match')
    parser.add_argument('--sample_size_perc', type=int, default=10)
    parser.add_argument('--sample_size_seed', type=int, default=0)
    parser.add_argument('--harm', action='store_true', default=False)

    args, _ = parser.parse_known_args()
    return args


def goalDNN_unharmed_input_gen(input_path,
                               output_path,
                               dataset_pair,
                               train_names=['unmatch2match_train'],
                               val_names=['unmatch2match_val'],
                               test_names=['unmatch2match_test'],
                               nb_folds=10):
    """
    Generate input for goalDNN using unharmonized data

    Args:
        input_path (str): Path for input data
        output_path (str): Path for output
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        train_names (list, optional): Names for training data.
        val_names (list, optional): Names for val data.
        test_names (list, optional): Names for test data
        nb_folds (int, optional): _description_. Defaults to 10.
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    numerical_features = ROIs + ['AGE', 'SITE', 'SEX', 'MMSE', 'DX']
    dataset = (dataset_pair.split('-'))[1]
    for fold in range(nb_folds):
        fold_input_path = os.path.join(input_path, str(fold))
        fold_output_path = os.path.join(output_path, str(fold))
        create_folder(fold_output_path)
        # load data
        train_file_names = [train_name + '.csv' for train_name in train_names]
        train_df = csvs2df(fold_input_path, train_file_names)
        train_df = train_df[cols]
        val_file_names = [val_name + '.csv' for val_name in val_names]
        val_df = csvs2df(fold_input_path, val_file_names)
        val_df = val_df[cols]
        test_file_names = [test_name + '.csv' for test_name in test_names]
        test_df = csvs2df(fold_input_path, test_file_names)
        test_df = test_df[cols]
        # extract site-wise data
        train_data = train_df[train_df.SITE == 0]  # unmatched ADNI for train
        val_data = val_df[val_df.SITE == 0]  # unmatched ADNI for tunning
        test_data_ADNI = test_df[test_df.SITE == 0]  # matched ADNI for testing
        test_data_dataset = test_df[test_df.SITE == 1]  # matched AIBL/MACC
        # processing for train_data
        train = dict()
        train_data.loc[:, 'SEX'] -= 1.0
        train['mean'] = train_data[numerical_features].mean()
        train['std'] = train_data[numerical_features].std()
        # forward filling
        train_data = forward_filling(train_data, 'MMSE')
        train_data = forward_filling(train_data, 'DX')
        # fill 0 for missing ROIs
        train_data.fillna(value=0, inplace=True)
        train['mean'][['SEX', 'DX', 'SITE']] = 0.
        train['std'][['SEX', 'DX', 'SITE']] = 1.
        train['data'] = z_norm(train_data[numerical_features], train['mean'],
                               train['std'])
        train['RID'] = train_data[['RID']]
        save_pkl(train, os.path.join(fold_output_path, 'train.pkl'))
        # preprocessing for val_data
        val = dict()
        val['mean'] = train['mean']
        val['std'] = train['std']
        # we need to drop tps without MMSE or DX
        val_data.dropna(subset=['DX', 'MMSE'], inplace=True)
        val_data.fillna(value=0, inplace=True)
        val['data'] = z_norm(val_data[numerical_features], val['mean'],
                             val['std'])
        val['RID'] = val_data[['RID']]
        save_pkl(val, os.path.join(fold_output_path, 'val.pkl'))
        # processing for test_data_ADNI
        test_ADNI = dict()
        test_ADNI['mean'] = train['mean']
        test_ADNI['std'] = train['std']
        test_data_ADNI.dropna(subset=['DX', 'MMSE'], inplace=True)
        test_data_ADNI.fillna(value=0, inplace=True)
        test_ADNI['data'] = z_norm(test_data_ADNI[numerical_features],
                                   test_ADNI['mean'], test_ADNI['std'])
        test_ADNI['RID'] = test_data_ADNI[['RID']]
        save_pkl(test_ADNI,
                 os.path.join(fold_output_path, 'test_ADNI_unharm.pkl'))
        # processing for test_data_dataset
        test_dataset = dict()
        test_dataset['mean'] = train['mean']
        test_dataset['std'] = train['std']
        test_data_dataset.dropna(subset=['DX', 'MMSE'], inplace=True)
        test_data_dataset.fillna(value=0, inplace=True)
        test_dataset['data'] = z_norm(test_data_dataset[numerical_features],
                                      test_dataset['mean'],
                                      test_dataset['std'])
        test_dataset['RID'] = test_data_dataset[['RID']]
        save_pkl(
            test_dataset,
            os.path.join(fold_output_path, 'test_' + dataset + '_unharm.pkl'))


def goalDNN_harmed_input_gen(input_path,
                             output_path,
                             dataset_pair,
                             test_name,
                             harm_model,
                             nb_folds=10):
    """
    Generate input for goalDNN using harmonized data by <harm_model>

    Args:
        input_path (str): Path for input data
        output_path (str): Path for output
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        test_name (str): Name for test data
        harm_model (str): Name for harmonziation model
        nb_folds (int, optional): Number of folds. Defaults to 10.
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    numerical_features = ROIs + ['AGE', 'SITE', 'SEX', 'MMSE', 'DX']
    dataset = (dataset_pair.split('-'))[1]
    for fold in range(nb_folds):
        fold_input_path = os.path.join(input_path, str(fold))
        fold_output_path = os.path.join(output_path, str(fold))
        create_folder(fold_output_path)
        # load data
        train = load_pkl(os.path.join(fold_output_path, 'train.pkl'))
        mean = train['mean']
        std = train['std']
        test_df = pd.read_csv(
            os.path.join(fold_input_path, test_name + '.csv'))
        test_df = test_df[cols]
        test_data_dataset = test_df[test_df.SITE == 1]
        test_dataset = dict()
        test_dataset['mean'] = mean
        test_dataset['std'] = std
        test_data_dataset.dropna(subset=['DX', 'MMSE'], inplace=True)
        test_data_dataset.fillna(value=0, inplace=True)
        test_dataset['data'] = z_norm(test_data_dataset[numerical_features],
                                      test_dataset['mean'],
                                      test_dataset['std'])
        test_dataset['RID'] = test_data_dataset[['RID']]
        save_pkl(
            test_dataset,
            os.path.join(fold_output_path,
                         'test_' + dataset + harm_model + '_harm.pkl'))


def VAE_input_gen(input_path,
                  output_path,
                  train_name='unmatch2match_train',
                  val_name='unmatch2match_val',
                  test_names=['unmatch2match_test'],
                  nb_folds=10):
    """
    Generate input for cVAE and gcVAE models

    Args:
        input_path (str): Path for input data
        output_path (str): Path for output
        train_name (str, optional): Name of training csv. Defaults to 'train'.
        val_name (str, optional): Name of val csv. Defaults to 'val'.
        test_names (list, optional): List of name of test csv.
        nb_folds (int, optional): Number of folds. Defaults to 10.
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    numerical_features = ROIs + ['AGE', 'SITE', 'SEX', 'MMSE', 'DX']
    for fold in range(nb_folds):
        fold_input_path = os.path.join(input_path, str(fold))
        fold_output_path = os.path.join(output_path, str(fold))
        create_folder(fold_output_path)
        # load data
        train_file_name = train_name + '.csv'
        train_df = pd.read_csv(os.path.join(fold_input_path, train_file_name))
        train_df = train_df[cols]
        val_file_name = val_name + '.csv'
        val_df = pd.read_csv(os.path.join(fold_input_path, val_file_name))
        val_df = val_df[cols]

        # extract for train val test
        train_data = train_df  # unmatched for train
        val_data = copy.deepcopy(val_df)  # unmatched for tunning

        # processing for train_data
        train = dict()
        train_data.loc[:, 'SEX'] -= 1.0
        train['mean'] = train_data[numerical_features].mean()
        train['std'] = train_data[numerical_features].std()
        # note this is for TaskcVAE part
        # even cVAE is not using MMSE and DX, taskcVAE would use them to
        # train model !!!
        train_data = forward_filling(train_data, 'MMSE')
        train_data = forward_filling(train_data, 'DX')
        # fill 0 for missing ROIs
        train_data.fillna(value=0, inplace=True)
        train['mean'][['SEX', 'DX', 'SITE']] = 0.
        train['std'][['SEX', 'DX', 'SITE']] = 1.
        train['data'] = z_norm(train_data[numerical_features], train['mean'],
                               train['std'])
        train['RID'] = train_data[['RID']]
        save_pkl(train, os.path.join(fold_output_path, 'train.pkl'))
        # preprocessing for val_data
        val = dict()
        val['mean'] = train['mean']
        val['std'] = train['std']
        val_data.fillna(value=0, inplace=True)
        val['data'] = z_norm(val_data[numerical_features], val['mean'],
                             val['std'])
        val['RID'] = val_data[['RID']]
        save_pkl(val, os.path.join(fold_output_path, 'val.pkl'))
        val = dict()
        val['mean'] = train['mean']
        val['std'] = train['std']
        val_data = copy.deepcopy(val_df)
        val_data.dropna(subset=['DX', 'MMSE'], inplace=True)
        val_data.fillna(value=0, inplace=True)
        val['data'] = z_norm(val_data[numerical_features], val['mean'],
                             val['std'])
        val['RID'] = val_data[['RID']]
        save_pkl(val, os.path.join(fold_output_path, 'val_gcVAE.pkl'))
        # in case of multiple test files
        for test_name in test_names:
            test_file_name = test_name + '.csv'
            test_df = pd.read_csv(
                os.path.join(fold_input_path, test_file_name))
            test_df = test_df[cols]
            # processing for test_data
            test_data = test_df  # matched for testing
            test = dict()
            test['mean'] = train['mean']
            test['std'] = train['std']
            test_data.fillna(value=0, inplace=True)
            test['data'] = z_norm(test_data[numerical_features], test['mean'],
                                  test['std'])
            test['RID'] = test_data[['RID']]
            save_pkl(test, os.path.join(fold_output_path, test_name + '.pkl'))


def goalDNN_harmed_input_gen_wrapper(goalDNN_input_path, harm_output_path,
                                     dataset_pair, test_name):
    """
    Generate ComBat/cVAE/gcVAE harmed input for goalDNN

    Args:
        goalDNN_input_path (str): Path for goalDNN input data
        harm_output_path (str): Path for harmonizaion output
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
    """
    # cVAE
    model_harm_output_path = os.path.join(harm_output_path, 'cVAE')
    goalDNN_harmed_input_gen(model_harm_output_path, goalDNN_input_path,
                             dataset_pair, test_name, 'cVAE')
    # gcVAE
    model_harm_output_path = os.path.join(harm_output_path, 'gcVAE')
    goalDNN_harmed_input_gen(model_harm_output_path, goalDNN_input_path,
                             dataset_pair, test_name, 'gcVAE')
    combat_test_name = test_name.split('-')[0]
    # ComBat
    model_harm_output_path = os.path.join(harm_output_path, 'ComBat',
                                          'AGE_SEX')
    goalDNN_harmed_input_gen(model_harm_output_path, goalDNN_input_path,
                             dataset_pair, combat_test_name, 'ComBat')
    # ComBat4cov
    model_harm_output_path = os.path.join(harm_output_path, 'ComBat',
                                          'AGE_SEX_MMSE_DX')
    goalDNN_harmed_input_gen(model_harm_output_path, goalDNN_input_path,
                             dataset_pair, combat_test_name, 'ComBat4cov')


def goalDNN_harmed_input_gen_top_wrapper(data_path,
                                         dataset_pair,
                                         exp,
                                         sample_size_perc=10,
                                         sample_size_seed=0):
    """
    Top wrapper function for generating harmed input for goalDNN

    Args:
        data_path (str): Path for data
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        exp (str): Name of experiments
        sample_size_perc (int, optional): Percentage for sampling.
        sample_size_seed (int, optional): Random seed for samping.
    """
    if exp == 'sample_size':
        goalDNN_input_path = os.path.join(data_path, 'goalDNN_input',
                                          dataset_pair)
        harm_output_path = os.path.join(data_path, 'harm_output', dataset_pair)
        test_name = 'unmatch2match_test-map2ADNI'
        goalDNN_harmed_input_gen_wrapper(goalDNN_input_path, harm_output_path,
                                         dataset_pair, test_name)
    else:
        exp_out_path = os.path.join(data_path, exp)
        goalDNN_input_path = os.path.join(exp_out_path, 'goalDNN_input',
                                          dataset_pair)
        harm_output_path = os.path.join(exp_out_path, 'harm_output',
                                        dataset_pair)
        test_name = exp + '_test-map2ADNI'
        goalDNN_harmed_input_gen_wrapper(goalDNN_input_path, harm_output_path,
                                         dataset_pair, test_name)


def goalDNN_unharmed_input_gen_top_wrapper(data_path,
                                           dataset_pair,
                                           exp,
                                           sample_size_perc=10,
                                           sample_size_seed=0):
    """
    Top wrapper function for generating unharmed input for goalDNN

    Args:
        data_path (str): Path for data
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        exp (str): Name of experiments
        sample_size_perc (int, optional): Percentage for sampling.
        sample_size_seed (int, optional): Random seed for samping.
    """
    exp_out_path = os.path.join(data_path, exp)
    if exp == 'sample_size':
        out_folder = str(sample_size_perc) + 'perc_seed' + \
            str(sample_size_seed)
        goalDNN_input_path = os.path.join(exp_out_path, out_folder,
                                          'goalDNN_input', dataset_pair)
        splits_path = os.path.join(exp_out_path, out_folder, 'splits',
                                   dataset_pair)
        goalDNN_unharmed_input_gen(
            splits_path,
            goalDNN_input_path,
            dataset_pair,
            train_names=['unmatch2match_train'],
            val_names=['unmatch2match_val'],
            test_names=['unmatch2match_test'])
    elif exp == 'unmatch2match':
        goalDNN_input_path = os.path.join(exp_out_path, 'goalDNN_input',
                                          dataset_pair)
        splits_path = os.path.join(data_path, 'splits', dataset_pair)
        goalDNN_unharmed_input_gen(
            splits_path,
            goalDNN_input_path,
            dataset_pair,
            train_names=[exp + '_train'],
            val_names=[exp + '_val'],
            test_names=[exp + '_test'])
    else:
        goalDNN_input_path = os.path.join(exp_out_path, 'goalDNN_input',
                                          dataset_pair)
        splits_path = os.path.join(data_path, 'splits', dataset_pair)
        goalDNN_unharmed_input_gen(
            splits_path,
            goalDNN_input_path,
            dataset_pair,
            train_names=['match2unmatch_train', 'unmatch2match_train'],
            val_names=['match2unmatch_val', 'unmatch2match_val'],
            test_names=[exp + '_test'])


def VAE_input_gen_top_wrapper(data_path,
                              dataset_pair,
                              exp,
                              sample_size_perc=10,
                              sample_size_seed=0):
    """
    Top wrapper function for generating input for VAE models

    Args:
        data_path (str): Path for data
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        exp (str): Name of experiments
        sample_size_perc (int, optional): Percentage for sampling.
        sample_size_seed (int, optional): Random seed for samping.
    """
    exp_out_path = os.path.join(data_path, exp)
    if exp == 'sample_size':
        out_folder = str(sample_size_perc) + 'perc_seed' + \
            str(sample_size_seed)
        test_names = [
            'unmatch2match_train_full', 'unmatch2match_val_full',
            'unmatch2match_test'
        ]
        harm_input_path = os.path.join(exp_out_path, out_folder, 'harm_input',
                                       dataset_pair)
        splits_path = os.path.join(exp_out_path, out_folder, 'splits',
                                   dataset_pair)
        VAE_input_gen(
            splits_path,
            harm_input_path,
            train_name='unmatch2match_train',
            val_name='unmatch2match_val',
            test_names=test_names)
    else:
        harm_input_path = os.path.join(exp_out_path, 'harm_input',
                                       dataset_pair)
        splits_path = os.path.join(data_path, 'splits', dataset_pair)
        VAE_input_gen(
            splits_path,
            harm_input_path,
            train_name=exp + '_train',
            val_name=exp + '_val',
            test_names=[exp + '_test'])


def wrapper(args):
    """
    Wrapper function for generating input for An2022_gcVAE

    Args:
        args (tuple): Parameters
    """
    assert args.exp in ['match2unmatch', 'unmatch2match', 'sample_size']
    assert args.dataset_pair in ['ADNI-AIBL', 'ADNI-MACC']
    if args.harm:
        goalDNN_harmed_input_gen_top_wrapper(args.data_path, args.dataset_pair,
                                             args.exp, args.sample_size_perc,
                                             args.sample_size_seed)
    else:
        # generate unharmonized input for VAE and goalDNN
        goalDNN_unharmed_input_gen_top_wrapper(
            args.data_path, args.dataset_pair, args.exp, args.sample_size_perc,
            args.sample_size_seed)
        VAE_input_gen_top_wrapper(args.data_path, args.dataset_pair, args.exp,
                                  args.sample_size_perc, args.sample_size_seed)


if __name__ == '__main__':
    wrapper(input_gen_args_parser())
