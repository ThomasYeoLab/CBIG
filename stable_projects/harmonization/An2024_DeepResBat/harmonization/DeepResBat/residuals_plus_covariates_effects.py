#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import pandas as pd
from typing import List
from config import global_config
from utils.misc import txt2list


def harmonized_res_plus_covar_effects_args_parser():
    parser = argparse.ArgumentParser(prog='G+FArgs')
    parser.add_argument('--harm_output_path', type=str, default='/')
    parser.add_argument('--nb_folds', type=int, default=10)
    parser.add_argument('--prefix', type=str, default='unmatch2match_')
    parser.add_argument('--g_files',
                        type=List,
                        default=['-recon', '-intermediate', '-map2ADNI'])
    parser.add_argument('--sets', type=List, default=['train', 'val', 'test'])
    parser.add_argument('--f_sufix', type=str, default='_F')
    parser.add_argument('--g_sufix', type=str, default='_G')
    parser.add_argument('--model', type=str, default='TrueBat')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')

    args, _ = parser.parse_known_args()
    return args


def harmonized_res_plus_covar_effects(args):
    ROIs = txt2list(global_config.ROI_features_path)
    cols = txt2list(global_config.columns_path)
    for fold in range(args.nb_folds):
        fold_output_path = os.path.join(args.harm_output_path, args.model,
                                        args.dataset_pair, str(fold))
        for s in args.sets:
            # load covariates effects file
            f_df = pd.read_csv(
                os.path.join(fold_output_path,
                             args.prefix + s + args.f_sufix + '.csv'))
            f_df = f_df[cols]
            for g_file in args.g_files:
                g_df = pd.read_csv(
                    os.path.join(
                        fold_output_path,
                        args.prefix + s + g_file + args.g_sufix + '.csv'))
                g_df = g_df[cols]
                g_df[ROIs] = g_df[ROIs] + f_df[ROIs]
                save_path = os.path.join(fold_output_path,
                                         args.prefix + s + g_file + '.csv')
                g_df.to_csv(save_path, index=False, sep=',')


if __name__ == '__main__':
    harmonized_res_plus_covar_effects(
        harmonized_res_plus_covar_effects_args_parser())
