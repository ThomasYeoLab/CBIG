#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import shutil
import itertools
import pandas as pd
from utils.misc import create_nfolds_folder, create_folder


def sample_data_args_parser():
    """
    Parameters for sampling data
    """
    parser = argparse.ArgumentParser(prog='SampleDataArgs')
    parser.add_argument('--root_path', type=str, default='/')
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')

    args, _ = parser.parse_known_args()
    return args


def prepare_folders_wrapper(root_path, dataset_pair, nb_seeds=10):
    """
    Create folders for 9 * 10 experiments with different sample size

    Args:
        root_path (str): Root path for An2022_gcVAE
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        nb_seeds (int, optional): Number of random seeds. Defaults to 10.
    """
    portions = [str(int(perc * 10)) + 'perc' for perc in range(1, 10)]
    seeds = ['_seed' + str(seed) for seed in range(10, nb_seeds + 10)]
    out_folders = [
        "".join(elem) for elem in itertools.product(portions, seeds)
    ]
    data_path = os.path.join(root_path, 'data', 'sample_size')
    checkpoints_path = os.path.join(root_path, 'checkpoints', 'sample_size')
    results_path = os.path.join(root_path, 'results', 'sample_size')
    for out_folder in out_folders:
        prepare_folders(data_path, checkpoints_path, results_path, out_folder,
                        dataset_pair)


def prepare_folders(data_path, checkpoints_path, results_path, out_folder,
                    dataset_pair):
    """
    Create folders for single <out_folder>

    Args:
        data_path (str): Path for data
        checkpoints_path (str): Path for checkpoints
        results_path (str): Path for results
        out_folder (str): Path for output
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
    """
    # data
    create_nfolds_folder(
        os.path.join(data_path, out_folder, 'splits', dataset_pair))
    create_nfolds_folder(
        os.path.join(data_path, out_folder, 'harm_input', dataset_pair))
    create_nfolds_folder(
        os.path.join(data_path, out_folder, 'goalDNN_input', dataset_pair))

    create_folder(
        os.path.join(data_path, out_folder, 'harm_output', dataset_pair,
                     'ComBat'))
    create_nfolds_folder(
        os.path.join(data_path, out_folder, 'harm_output', dataset_pair,
                     'cVAE'))
    create_nfolds_folder(
        os.path.join(data_path, out_folder, 'harm_output', dataset_pair,
                     'gcVAE'))

    create_folder(
        os.path.join(checkpoints_path, out_folder, 'eval_model', 'goalDNN',
                     dataset_pair))

    create_folder(
        os.path.join(checkpoints_path, out_folder, 'eval_model', 'XGBoost',
                     dataset_pair))
    create_folder(
        os.path.join(checkpoints_path, out_folder, 'harm_model', 'cVAE',
                     dataset_pair))

    create_folder(
        os.path.join(checkpoints_path, out_folder, 'harm_model', 'gcVAE',
                     dataset_pair))


def sample_portion_data(data_path, dataset_pair, nb_seeds=10, nb_folds=10):
    """
    Sample portion training + val data

    Args:
        data_path (str): Pathfor data
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
        nb_seeds (int, optional): Number of random seeds. Defaults to 10.
        nb_folds (int, optional): Number of random seeds. Defaults to 10.
    """
    percs = [90, 80, 70, 60, 50, 40, 30, 20, 10]
    origin_splits_path = os.path.join(data_path, 'splits', dataset_pair)
    origin_train_csv = 'unmatch2match_train.csv'
    origin_val_csv = 'unmatch2match_val.csv'
    origin_test_csv = 'unmatch2match_test.csv'
    output_path = os.path.join(data_path, 'sample_size')
    for fold in range(nb_folds):
        fold_origin_splits_path = os.path.join(origin_splits_path, str(fold))
        for seed in range(nb_seeds, 10 + nb_seeds):
            # read data
            train_df = pd.read_csv(
                os.path.join(fold_origin_splits_path, origin_train_csv))
            val_df = pd.read_csv(
                os.path.join(fold_origin_splits_path, origin_val_csv))
            train_df_adni = train_df[train_df.SITE == 0]
            train_df_nonadni = train_df[train_df.SITE == 1]
            val_df_adni = val_df[val_df.SITE == 0]
            val_df_nonadni = val_df[val_df.SITE == 1]
            init_frac_ratio = percs[0] / 100
            for i, perc in enumerate(percs):
                frac_ratio = perc / 100
                out_folder = str(perc) + 'perc_seed' + str(seed)
                frac_output_path = os.path.join(output_path, out_folder,
                                                'splits', dataset_pair)
                if i == 0:
                    # firstly sampled 90%
                    # train
                    sampled_train_df_adni_90 = train_df_adni.sample(
                        frac=init_frac_ratio, replace=False, random_state=seed)
                    sampled_train_df_nonadni_90 = train_df_nonadni.sample(
                        frac=init_frac_ratio, replace=False, random_state=seed)
                    sampled_train_df_90 = sampled_train_df_adni_90.append(
                        sampled_train_df_nonadni_90, ignore_index=True)
                    # save train
                    sampled_train_df_90.to_csv(
                        os.path.join(frac_output_path, str(fold),
                                     origin_train_csv),
                        index=False,
                        sep=',')
                    # val
                    sampled_val_df_adni_90 = val_df_adni.sample(
                        frac=init_frac_ratio, replace=False, random_state=seed)
                    sampled_val_df_nonadni_90 = val_df_nonadni.sample(
                        frac=init_frac_ratio, replace=False, random_state=seed)
                    sampled_val_df_90 = sampled_val_df_adni_90.append(
                        sampled_val_df_nonadni_90, ignore_index=True)
                    # save val
                    sampled_val_df_90.to_csv(
                        os.path.join(frac_output_path, str(fold),
                                     origin_val_csv),
                        index=False,
                        sep=',')
                    # copy test data & full training + val
                    shutil.copyfile(
                        src=os.path.join(fold_origin_splits_path,
                                         origin_test_csv),
                        dst=os.path.join(frac_output_path, str(fold),
                                         origin_test_csv))
                    shutil.copyfile(
                        src=os.path.join(fold_origin_splits_path,
                                         origin_train_csv),
                        dst=os.path.join(frac_output_path, str(fold),
                                         'unmatch2match_train_full.csv'))
                    shutil.copyfile(
                        src=os.path.join(fold_origin_splits_path,
                                         origin_val_csv),
                        dst=os.path.join(frac_output_path, str(fold),
                                         'unmatch2match_val_full.csv'))
                else:
                    # sampled
                    sampled_train_df_adni = sampled_train_df_adni_90.sample(
                        frac=frac_ratio / init_frac_ratio,
                        replace=False,
                        random_state=seed)
                    sampled_train_df_nonadni = \
                        sampled_train_df_nonadni_90.sample(
                            frac=frac_ratio / init_frac_ratio,
                            replace=False,
                            random_state=seed)
                    sampled_train_df = sampled_train_df_adni.append(
                        sampled_train_df_nonadni, ignore_index=True)
                    sampled_train_df_adni_90 = sampled_train_df_adni
                    sampled_train_df_nonadni_90 = sampled_train_df_nonadni

                    # save train
                    sampled_train_df.to_csv(
                        os.path.join(frac_output_path, str(fold),
                                     origin_train_csv),
                        index=False,
                        sep=',')
                    # val
                    sampled_val_df_adni = sampled_val_df_adni_90.sample(
                        frac=frac_ratio / init_frac_ratio,
                        replace=False,
                        random_state=seed)
                    sampled_val_df_nonadni = sampled_val_df_nonadni_90.sample(
                        frac=frac_ratio / init_frac_ratio,
                        replace=False,
                        random_state=seed)
                    sampled_val_df = sampled_val_df_adni.append(
                        sampled_val_df_nonadni, ignore_index=True)
                    # update
                    sampled_val_df_adni_90 = sampled_val_df_adni
                    sampled_val_df_nonadni_90 = sampled_val_df_nonadni
                    init_frac_ratio = frac_ratio
                    # save val
                    sampled_val_df.to_csv(
                        os.path.join(frac_output_path, str(fold),
                                     origin_val_csv),
                        index=False,
                        sep=',')
                    # copy test data & full training + val
                    shutil.copyfile(
                        src=os.path.join(fold_origin_splits_path,
                                         origin_test_csv),
                        dst=os.path.join(frac_output_path, str(fold),
                                         origin_test_csv))
                    shutil.copyfile(
                        src=os.path.join(fold_origin_splits_path,
                                         origin_train_csv),
                        dst=os.path.join(frac_output_path, str(fold),
                                         'unmatch2match_train_full.csv'))
                    shutil.copyfile(
                        src=os.path.join(fold_origin_splits_path,
                                         origin_val_csv),
                        dst=os.path.join(frac_output_path, str(fold),
                                         'unmatch2match_val_full.csv'))


def copy_optimal_hyperparameters(root_path, dataset_pair):
    """
    Copy optimial hyper-parameters configureations to <out_folder>

    Args:
        root_path (str): Root path for An2022_gcVAE
        dataset_pair (str): ADNI-AIBL or ADNI-MACC
    """
    checkpoints_path = os.path.join(root_path, 'checkpoints')
    unmatch2match_checkpoints_path = os.path.join(checkpoints_path,
                                                  'unmatch2match')
    sample_size_checkpoints_path = os.path.join(checkpoints_path,
                                                'sample_size')

    portions = [str(int(perc * 10)) + 'perc' for perc in range(1, 10)]
    seeds = ['_seed' + str(seed) for seed in range(10, 20)]
    out_folders = [
        "".join(elem) for elem in itertools.product(portions, seeds)
    ]
    for out_folder in out_folders:
        out_folder_checkpoints_path = os.path.join(
            sample_size_checkpoints_path, out_folder)
        # goalDNN
        shutil.copy(
            src=os.path.join(unmatch2match_checkpoints_path, 'eval_model',
                             'goalDNN', dataset_pair, 'hord_params.csv'),
            dst=os.path.join(out_folder_checkpoints_path, 'eval_model',
                             'goalDNN', dataset_pair, 'hord_params.csv'))
        # cVAE
        shutil.copy(
            src=os.path.join(unmatch2match_checkpoints_path, 'harm_model',
                             'cVAE', dataset_pair, 'hord_params.csv'),
            dst=os.path.join(out_folder_checkpoints_path, 'harm_model', 'cVAE',
                             dataset_pair, 'hord_params.csv'))
        # gcVAE
        shutil.copy(
            src=os.path.join(unmatch2match_checkpoints_path, 'harm_model',
                             'gcVAE', dataset_pair, 'grid_params.csv'),
            dst=os.path.join(out_folder_checkpoints_path, 'harm_model',
                             'gcVAE', dataset_pair, 'grid_params.csv'))


def wrapper(args):
    """
    Wrapper function for sampling data

    Args:
        args (tuple): Parameters
    """
    prepare_folders_wrapper(args.root_path, args.dataset_pair)
    copy_optimal_hyperparameters(args.root_path, args.dataset_pair)
    sample_portion_data(args.data_path, args.dataset_pair)


if __name__ == '__main__':
    wrapper(sample_data_args_parser())
