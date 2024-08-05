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
from utils.input_generation import deepharm_input_gen
from utils.misc import txt2list, save_df, clean_nfolds_files


def residual_gen_args_parser():
    """
    Feeding parameters
    """
    parser = argparse.ArgumentParser(prog='ResidualGeneratorArgs')
    parser.add_argument('--dataset_pair', type=str, default='/')
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--harm_input_path', type=str, default='/')
    parser.add_argument('--harm_output_path', type=str, default='/')
    parser.add_argument('--model', type=str, default='DeepResBat')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--f_sufix', type=str, default='_F')
    parser.add_argument('--ysf_sufix', type=str, default='_G')
    parser.add_argument('--train_name',
                        type=str,
                        default='unmatch2match_train')
    parser.add_argument('--val_name', type=str, default='unmatch2match_val')
    parser.add_argument('--test_name', type=str, default='unmatch2match_test')
    parser.add_argument('--demean', action='store_false', default=False)

    args, _ = parser.parse_known_args()
    return args


def generate_residual(org_csv_path, covs_csvs_path, cols, ROIs, res_csv_path):
    """
    Generate residuals
    """
    org_df = pd.read_csv(org_csv_path)
    org_df = org_df[cols]
    covs_df = pd.read_csv(covs_csvs_path)
    covs_df = covs_df[cols]
    # substract
    res_df = copy.deepcopy(org_df)
    res_df[ROIs] = org_df[ROIs] - covs_df[ROIs]
    save_df(res_df, res_csv_path)


def residual_gen_wrapper(args):
    """
    Wrapper function for generating residuals
    """
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    # 1. get Y - F(X)
    for fold in range(args.nb_folds):
        fold_data_path = os.path.join(args.data_path, args.dataset_pair,
                                      str(fold))
        fold_harm_output_path = os.path.join(args.harm_output_path, args.model,
                                             args.dataset_pair, str(fold))
        # train
        generate_residual(
            os.path.join(fold_data_path, args.train_name + '.csv'),
            os.path.join(fold_harm_output_path,
                         args.train_name + args.f_sufix + '.csv'), cols, ROIs,
            os.path.join(fold_harm_output_path,
                         args.train_name + args.ysf_sufix + '.csv'))
        # val
        generate_residual(
            os.path.join(fold_data_path, args.val_name + '.csv'),
            os.path.join(fold_harm_output_path,
                         args.val_name + args.f_sufix + '.csv'), cols, ROIs,
            os.path.join(fold_harm_output_path,
                         args.val_name + args.ysf_sufix + '.csv'))
        # test
        generate_residual(
            os.path.join(fold_data_path, args.test_name + '.csv'),
            os.path.join(fold_harm_output_path,
                         args.test_name + args.f_sufix + '.csv'), cols, ROIs,
            os.path.join(fold_harm_output_path,
                         args.test_name + args.ysf_sufix + '.csv'))
    # 2. call vae_input_gen to get residual input for deep learning models
    deepharm_input_args = copy.deepcopy(args)
    deepharm_input_args.train_name = args.train_name + args.ysf_sufix
    deepharm_input_args.val_name = args.val_name + args.ysf_sufix
    deepharm_input_args.test_name = args.test_name + args.ysf_sufix
    deepharm_input_args.input_path = os.path.join(args.harm_output_path,
                                                  args.model)
    deepharm_input_gen(deepharm_input_args)
    # 3. clean the intermediate files
    clean_nfolds_files(
        os.path.join(args.harm_output_path, args.model, args.dataset_pair), [
            args.train_name + args.ysf_sufix + '.csv',
            args.val_name + args.ysf_sufix + '.csv',
            args.test_name + args.ysf_sufix + '.csv',
        ])


if __name__ == '__main__':
    residual_gen_wrapper(residual_gen_args_parser())
