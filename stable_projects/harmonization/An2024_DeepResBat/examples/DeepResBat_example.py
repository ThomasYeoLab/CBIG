#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Example code for applying DeepResBat harmonization
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import argparse
import numpy as np
import pandas as pd
from sklearn.datasets import make_regression
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from utils.misc import create_folder, txt2list, save_df
from utils.input_generation import deepharm_input_gen, \
    deepharm_input_gen_args_parser
from config import global_config

from harmonization.DeepResBat.covariates_effects_estimator import \
    covariates_effects_estimator_parser, covariatess_effects_estimator
from harmonization.DeepResBat.residuals_generator import \
    residual_gen_args_parser, residual_gen_wrapper
from harmonization.DeepResBat.residuals_harmonizer_fit import \
    residuals_harmonizer_fit_args_parser, residuals_harmonizer_fit
from harmonization.DeepResBat.residuals_harmonizer_infer import \
    residuals_harmonizer_infer_args_parser, residuals_harmonizer_infer
from harmonization.DeepResBat.residuals_plus_covariates_effects import \
    harmonized_res_plus_covar_effects_args_parser, \
    harmonized_res_plus_covar_effects
from evaluation.association_analysis.manova import \
    manova_args_parser, manova_wrapper


def example_args_parser():
    """
    Parameters for running example code of An2024_DeepResBat project
    """
    parser = argparse.ArgumentParser(prog='ExampleArgs')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--nb_sites', type=int, default=2)
    parser.add_argument('--nb_features', type=int, default=108)
    parser.add_argument('--nb_informative', type=int, default=54)
    parser.add_argument('--nb_samples_per_site', type=int, default=1000)
    parser.add_argument('--train_ratio', type=float, default=0.8)
    parser.add_argument('--test_ratio', type=int, default=0.1)
    parser.add_argument('--working_dir',
                        type=str,
                        default=os.path.join(global_config.root_path,
                                             'examples'))
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')

    parser.add_argument('--gen_data', action='store_true', default=False)
    parser.add_argument('--covar_effects', action='store_true', default=False)
    parser.add_argument('--gen_res', action='store_true', default=False)
    parser.add_argument('--HORD', action='store_true', default=False)
    parser.add_argument('--harm_res', action='store_true', default=False)
    parser.add_argument('--add_back', action='store_true', default=False)
    parser.add_argument('--manova', action='store_true', default=False)
    parser.add_argument('--print_results', action='store_true', default=False)
    parser.add_argument('--unittest', action='store_true', default=False)

    args, _ = parser.parse_known_args()
    return args


def gen_example_data_step1(args):
    """
    Generate two heterogeneous datatases using sklearn
    """
    # generate data using sklearn.datasts.make_regression
    # simulaute brain ROI thickness/volumes by FreeSurfer
    x, y = make_regression(n_samples=args.nb_samples_per_site * args.nb_sites,
                           n_features=args.nb_features,
                           n_informative=args.nb_informative,
                           n_targets=1,
                           noise=0.1,
                           random_state=args.seed)
    y = np.reshape(y, (y.shape[0], -1))
    # scale y to [0, 30] to mimic MMSE range
    scaler = MinMaxScaler(feature_range=(0, 30)).fit(y)
    y = scaler.transform(y)
    y = y.astype(int)

    # unbalanced y distribution across sites
    # site0 has lower y, site1 has higher y
    low_index = np.where(y[:, 0] <= 10)
    low_x, low_y = x[low_index], y[low_index]
    mid_index = np.where((y[:, 0] > 10) & (y[:, 0] <= 20))
    mid_x, mid_y = x[mid_index], y[mid_index]
    high_index = np.where((y[:, 0] > 20) & (y[:, 0] <= 30))
    high_x, high_y = x[high_index], y[high_index]
    x_site0 = np.concatenate((low_x[:int(low_x.shape[0] * 0.7), :],
                              mid_x[:int(mid_x.shape[0] * 0.6), :],
                              high_x[:int(high_x.shape[0] * 0.2), :]),
                             axis=0)
    x_site1 = np.concatenate((low_x[int(low_x.shape[0] * 0.7):, :],
                              mid_x[int(mid_x.shape[0] * 0.6):, :],
                              high_x[int(high_x.shape[0] * 0.2):, :]),
                             axis=0)
    y_site0 = np.concatenate((low_y[:int(low_x.shape[0] * 0.7), :],
                              mid_y[:int(mid_x.shape[0] * 0.6), :],
                              high_y[:int(high_x.shape[0] * 0.2), :]),
                             axis=0)
    y_site1 = np.concatenate((low_y[int(low_x.shape[0] * 0.7):, :],
                              mid_y[int(mid_x.shape[0] * 0.6):, :],
                              high_y[int(high_x.shape[0] * 0.2):, :]),
                             axis=0)
    # add Gaussian noise on x for each site ==> random site effect
    rng = np.random.default_rng(seed=args.seed)
    noise_site0 = rng.normal(loc=0.3, scale=1, size=x_site0.shape)
    rng = np.random.default_rng(seed=args.seed)
    noise_site1 = rng.normal(loc=-0.3, scale=2, size=x_site1.shape)
    rng = np.random.default_rng(seed=args.seed)
    x_site0 = x_site0 + noise_site0
    x_site1 = x_site1 + noise_site1

    return x_site0, y_site0, x_site1, y_site1


def gen_example_data_step2(args, x_site0, y_site0, x_site1, y_site1):
    """
    Generate train/val/test split
    """
    # for site 0
    x_site0_train, x_site0_test, y_site0_train, y_site0_test =\
        train_test_split(
            x_site0, y_site0,
            test_size=args.test_ratio,
            random_state=args.seed)
    x_site0_train, x_site0_val, y_site0_train, y_site0_val =\
        train_test_split(
            x_site0_train,
            y_site0_train,
            test_size=(args.test_ratio / args.train_ratio),
            random_state=args.seed)
    # for site 1
    x_site1_train, x_site1_test, y_site1_train, y_site1_test =\
        train_test_split(
            x_site1, y_site1,
            test_size=args.test_ratio,
            random_state=args.seed)
    x_site1_train, x_site1_val, y_site1_train, y_site1_val =\
        train_test_split(
            x_site1_train,
            y_site1_train,
            test_size=(args.test_ratio / args.train_ratio),
            random_state=args.seed)

    # concatenate
    # train
    x_train = np.concatenate((x_site0_train, x_site1_train), axis=0)
    y_train = np.concatenate((y_site0_train, y_site1_train), axis=0)
    train_site = np.concatenate(
        (np.zeros((y_site0_train.shape[0], 1)),
            np.ones((y_site1_train.shape[0], 1))),
        axis=0)
    train_yx = np.concatenate((y_train, x_train), axis=1)
    # val
    x_val = np.concatenate((x_site0_val, x_site1_val), axis=0)
    y_val = np.concatenate((y_site0_val, y_site1_val), axis=0)
    val_site = np.concatenate(
        (np.zeros((y_site0_val.shape[0], 1)),
            np.ones((y_site1_val.shape[0], 1))),
        axis=0)
    val_yx = np.concatenate((y_val, x_val), axis=1)
    # test
    x_test = np.concatenate((x_site0_test, x_site1_test), axis=0)
    y_test = np.concatenate((y_site0_test, y_site1_test), axis=0)
    test_site = np.concatenate((np.zeros(
        (y_site0_test.shape[0], 1)), np.ones((y_site1_test.shape[0], 1))),
        axis=0)
    test_yx = np.concatenate((y_test, x_test), axis=1)

    return train_yx, train_site, val_yx, val_site, test_yx, test_site


def gen_example_data_step3(train_yx, train_site, val_yx, val_site, test_yx,
                           test_site):
    """
    Generate dataframe as format of real ADNI data
    """
    cols = txt2list(global_config.columns_path)
    n_train, n_val = train_site.shape[0], val_site.shape[0]

    def gen_example_data_dummy_df(yx, site, subRID_begin, cols):
        """
        Generate dummy dataframe
        """
        nb_subs = yx.shape[0]
        rids = np.array(
            [i for i in range(subRID_begin, subRID_begin + nb_subs)]).reshape(
                (-1, 1))
        dates = np.array(['2022-01-01' for _ in range(nb_subs)]).reshape(
            (-1, 1))
        rng = np.random.default_rng(seed=0)
        dxs = (rng.integers(low=0, high=3, size=nb_subs)).reshape((-1, 1))
        dxs = dxs.astype(int)
        rng = np.random.default_rng(seed=0)
        ages = rng.random((nb_subs, 1))
        rng = np.random.default_rng(seed=0)
        sexs = (rng.integers(low=0, high=2, size=nb_subs)).reshape((-1, 1))
        sexs = sexs.astype(int) + 1
        site = site.astype(int)
        array = np.concatenate((rids, dates, dxs, site, ages, sexs, yx),
                               axis=1)
        df = pd.DataFrame(data=array, columns=cols)
        return df

    # for train
    train_df = gen_example_data_dummy_df(train_yx, train_site, 0, cols)
    # for val
    val_df = gen_example_data_dummy_df(val_yx, val_site, n_train, cols)
    # for test
    test_df = gen_example_data_dummy_df(test_yx, test_site, n_train + n_val,
                                        cols)

    return train_df, val_df, test_df


def gen_example_data(args):
    """
    Generate toy example data
    """
    data_folder = os.path.join(args.working_dir, 'data')
    splits_path = os.path.join(data_folder, 'splits')
    harm_input_path = os.path.join(data_folder, 'unmatch2match', 'harm_input')
    harm_output_path = os.path.join(data_folder, 'unmatch2match',
                                    'harm_output')
    create_folder(os.path.join(splits_path, args.dataset_pair, '0'))
    create_folder(os.path.join(harm_input_path, args.dataset_pair, '0'))
    create_folder(os.path.join(harm_output_path, args.dataset_pair, '0'))

    # step1, use sklearn to generate two heterogeneous datatases
    x_site0, y_site0, x_site1, y_site1 = gen_example_data_step1(args)
    # step2, split train/val/test
    train_yx, s_train, val_yx, s_val, test_yx, s_test =\
        gen_example_data_step2(args, x_site0, y_site0, x_site1, y_site1)
    # step3, generate dataframe as format of real ADNI data
    train_df, val_df, test_df = gen_example_data_step3(train_yx, s_train,
                                                       val_yx, s_val, test_yx,
                                                       s_test)
    # step4, save as csv format
    save_df(
        train_df,
        os.path.join(splits_path, args.dataset_pair, '0',
                     'unmatch2match_train.csv'))
    save_df(
        val_df,
        os.path.join(splits_path, args.dataset_pair, '0',
                     'unmatch2match_val.csv'))
    save_df(
        test_df,
        os.path.join(splits_path, args.dataset_pair, '0',
                     'unmatch2match_test.csv'))
    # step5, call utils.input_generation to generate data
    deepharm_input_gen_args = deepharm_input_gen_args_parser()
    deepharm_input_gen_args.input_path = splits_path
    deepharm_input_gen_args.harm_input_path = harm_input_path
    deepharm_input_gen_args.nb_folds = 1
    deepharm_input_gen_args.dataset_pair = args.dataset_pair
    deepharm_input_gen(deepharm_input_gen_args)


def example_wrapper(args):
    """
    Wrapper function for examples/unittests
    """
    if args.unittest:
        args.working_dir = os.path.join(global_config.root_path, 'unit_tests')

    data_folder = os.path.join(args.working_dir, 'data')
    splits_path = os.path.join(data_folder, 'splits')
    harm_input_path = os.path.join(data_folder, 'unmatch2match', 'harm_input')
    harm_output_path = os.path.join(data_folder, 'unmatch2match',
                                    'harm_output')
    checkpoint_path = os.path.join(args.working_dir, 'checkpoints',
                                   'unmatch2match', 'harm_model')
    create_folder(data_folder)
    create_folder(splits_path)
    create_folder(harm_input_path)
    create_folder(harm_output_path)
    create_folder(checkpoint_path)

    if args.gen_data:
        # generate toy datasets
        gen_example_data(args)

    if args.covar_effects:
        # esitimate effects of covariates using XGBoost
        covariates_effects_estimator_args = \
            covariates_effects_estimator_parser()
        covariates_effects_estimator_args.data_path = splits_path
        covariates_effects_estimator_args.checkpoint_path = checkpoint_path
        covariates_effects_estimator_args.hyper_params_path = checkpoint_path
        covariates_effects_estimator_args.output_path = harm_output_path
        covariates_effects_estimator_args.dataset_pair = args.dataset_pair
        covariates_effects_estimator_args.model = 'DeepResBat'
        covariates_effects_estimator_args.nb_folds = 1
        covariates_effects_estimator_args.verbose = False
        covariatess_effects_estimator(covariates_effects_estimator_args)

    if args.gen_res:
        # generate covariate-free residuals
        residual_gen_args = residual_gen_args_parser()
        residual_gen_args.data_path = splits_path
        residual_gen_args.harm_input_path = harm_input_path
        residual_gen_args.harm_output_path = harm_output_path
        residual_gen_args.dataset_pair = args.dataset_pair
        residual_gen_args.model = 'DeepResBat'
        residual_gen_args.nb_folds = 1
        residual_gen_wrapper(residual_gen_args)

    if args.HORD:
        # HORD search should be only done on GPU!
        # HORD search is wrapped as a shell script
        # please refer to $ROOTDIR/HORD/DeepResBat_HORD_example.sh
        HORD_script_path = os.path.join(global_config.root_path, 'HORD',
                                        'CBIG_DeepResBat_HORD_example.sh')
        os.system(HORD_script_path)

    if args.harm_res:
        # Harmonize covariate-free residuals using cVAE
        # Fit models
        residuals_harmonizer_fit_args = residuals_harmonizer_fit_args_parser()
        residuals_harmonizer_fit_args.cpu = True
        residuals_harmonizer_fit_args.epochs = 5
        residuals_harmonizer_fit_args.data_path = os.path.join(
            harm_input_path, args.dataset_pair, '0')
        residuals_harmonizer_fit_args.checkpoint_path = os.path.join(
            checkpoint_path, 'DeepResBat', args.dataset_pair, '0')
        residuals_harmonizer_fit_args.isSaving = True
        residuals_harmonizer_fit_args.model_name = 'DeepResBat'
        residuals_harmonizer_fit_args.sufix = '_G'
        residuals_harmonizer_fit_args.lr = 0.01
        residuals_harmonizer_fit_args.drop_out = 0.5
        residuals_harmonizer_fit_args.alpha = 1
        residuals_harmonizer_fit_args.lambda_ = 1
        residuals_harmonizer_fit_args.gamma = 5
        residuals_harmonizer_fit_args.lr_step = 2
        residuals_harmonizer_fit_args.latent_dim = 64
        residuals_harmonizer_fit_args.nb_layers = 2
        residuals_harmonizer_fit_args.h1 = 256
        residuals_harmonizer_fit_args.h2 = 128
        residuals_harmonizer_fit(residuals_harmonizer_fit_args)
        # Infer
        parent_checkpoint_path = os.path.join(args.working_dir, 'checkpoints',
                                              'unmatch2match')
        residuals_harmonizer_infer_args = \
            residuals_harmonizer_infer_args_parser()
        residuals_harmonizer_infer_args.raw_data_path = splits_path
        residuals_harmonizer_infer_args.harm_input_path = harm_input_path
        residuals_harmonizer_infer_args.harm_output_path = harm_output_path
        residuals_harmonizer_infer_args.checkpoint_path = \
            parent_checkpoint_path
        residuals_harmonizer_infer_args.dataset_pair = args.dataset_pair
        residuals_harmonizer_infer_args.model_name = 'DeepResBat'
        residuals_harmonizer_infer_args.sufix = '_G'
        residuals_harmonizer_infer_args.nb_folds = 1
        residuals_harmonizer_infer(residuals_harmonizer_infer_args)
    if args.add_back:
        # add back estimated covariate effects to harmonized residuals
        harmonized_res_plus_covar_effects_args = \
            harmonized_res_plus_covar_effects_args_parser()
        harmonized_res_plus_covar_effects_args.harm_output_path = \
            harm_output_path
        harmonized_res_plus_covar_effects_args.dataset_pair = args.dataset_pair
        harmonized_res_plus_covar_effects_args.model = 'DeepResBat'
        harmonized_res_plus_covar_effects_args.g_sufix = '_G'
        harmonized_res_plus_covar_effects_args.nb_folds = 1
        harmonized_res_plus_covar_effects(
            harmonized_res_plus_covar_effects_args)
    if args.manova:
        manova_args = manova_args_parser()
        manova_args.input_path = os.path.join(harm_output_path, 'DeepResBat',
                                              args.dataset_pair)
        manova_results_path = os.path.join(args.working_dir, 'results',
                                           'assoc_manova')
        create_folder(manova_results_path)
        manova_args.output_path = manova_results_path
        manova_args.dataset_pair = args.dataset_pair
        manova_args.test_name = 'unmatch2match_test-intermediate'
        manova_args.model = 'DeepResBat'
        manova_args.type = 'MMSE'
        manova_args.nb_folds = 1
        manova_wrapper(manova_args)
    if args.print_results:
        manova_csv_path = os.path.join(manova_results_path, 'MMSE_demean',
                                       args.dataset_pair,
                                       'p_DeepResBat_ADNI-AIBL.csv')
        manova_df = pd.read_csv(manova_csv_path)
        manova_df.rename(columns={"P": "-log(P)"})
        print(manova_df)


if __name__ == '__main__':
    example_wrapper(example_args_parser())
