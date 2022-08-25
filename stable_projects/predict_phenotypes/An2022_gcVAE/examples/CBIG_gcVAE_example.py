#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
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
from utils.misc import create_folder, txt2list, save_df, myException
from config import global_config
from utils.input_generation\
    import goalDNN_unharmed_input_gen, goalDNN_harmed_input_gen, VAE_input_gen
from evaluation.goalDNN.train_goalDNN import train_goalDNN_args_parser
from evaluation.goalDNN.train_goalDNN import train as train_goaldnn
from evaluation.goalDNN.predict_goalDNN import predict_goalDNN_args_parser
from evaluation.goalDNN.predict_goalDNN import predict as predict_goaldnn
from evaluation.goalDNN.eval_prediction import eval_goalDNN_args_parser
from evaluation.goalDNN.eval_prediction import evaluation
from evaluation.XGBoost.train_XGBoost import train_XGBoost_args_parser
from evaluation.XGBoost.train_XGBoost import train as train_XGBoost
from harmonization.cVAE.train_cVAE import train_cvae_args_parser
from harmonization.cVAE.train_cVAE import train as train_cvae
from harmonization.cVAE.predict_cVAE import predict_cvae_args_parser
from harmonization.cVAE.predict_cVAE import predict as predict_cvae
from harmonization.gcVAE.train_gcVAE import train_gcVAE_args_parser
from harmonization.gcVAE.train_gcVAE import train as train_gcvae
from harmonization.gcVAE.predict_gcVAE import predict_gcvae_args_parser
from harmonization.gcVAE.predict_gcVAE import predict as predict_gcvae


def example_args_parser():
    """
    Parameters for running example code of An2022_gcVAE project
    """
    parser = argparse.ArgumentParser(prog='ExampleArgs')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--nb_sites', type=int, default=2)
    parser.add_argument('--nb_features', type=int, default=108)
    parser.add_argument('--nb_informative', type=int, default=54)
    parser.add_argument('--nb_samples_per_site', type=int, default=1000)
    parser.add_argument('--train_ratio', type=float, default=0.8)
    parser.add_argument('--test_ratio', type=int, default=0.1)
    parser.add_argument(
        '--working_dir',
        type=str,
        default=os.path.join(global_config.root_path, 'examples'))
    parser.add_argument('--model', type=str, default='goalDNN')
    parser.add_argument('--stage', type=str, default='prepare')
    parser.add_argument('--dataset_pair', type=str, default='ADNI-AIBL')
    parser.add_argument('--unittest', action='store_true', default=False)

    args, _ = parser.parse_known_args()
    return args


def gen_example_data_step1(args):
    """
    Generate two heterogeneous datatases using sklearn

    Args:
        args (tuple): Parameters
    """
    # generate data using sklearn.datasts.make_regression
    x, y = make_regression(
        n_samples=args.nb_samples_per_site * args.nb_sites,
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

    Args:
        args (tuple): Parmaters
        x_site0 (ndarray): X(features) from site 0
        y_site0 (ndarray): Y(traget) from site 0
        x_site1 (ndarray): X(features) from site 1
        y_site1 (ndarray): Y(traget) from site 0
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
    s_train = np.concatenate((np.zeros(
        (y_site0_train.shape[0], 1)), np.ones((y_site1_train.shape[0], 1))),
                             axis=0)
    train_yx = np.concatenate((y_train, x_train), axis=1)
    # val
    x_val = np.concatenate((x_site0_val, x_site1_val), axis=0)
    y_val = np.concatenate((y_site0_val, y_site1_val), axis=0)
    s_val = np.concatenate((np.zeros(
        (y_site0_val.shape[0], 1)), np.ones((y_site1_val.shape[0], 1))),
                           axis=0)
    val_yx = np.concatenate((y_val, x_val), axis=1)
    # test
    x_test = np.concatenate((x_site0_test, x_site1_test), axis=0)
    y_test = np.concatenate((y_site0_test, y_site1_test), axis=0)
    s_test = np.concatenate((np.zeros(
        (y_site0_test.shape[0], 1)), np.ones((y_site1_test.shape[0], 1))),
                            axis=0)
    test_yx = np.concatenate((y_test, x_test), axis=1)

    return train_yx, s_train, val_yx, s_val, test_yx, s_test


def gen_example_data_step3(train_yx, s_train, val_yx, s_val, test_yx, s_test):
    """
    generate dataframe as format of real ADNI data

    Args:
        train_yx (ndarray): [Y,X] for training
        s_train (ndarray): S(site vector) for training
        val_yx (ndarray): [Y,X] for validation
        s_val (ndarray): S(site vector) for validation
        test_yx (ndarray): [Y,X] for testing
        s_test (ndarray): S(site vector) for testing
    """
    cols = txt2list(global_config.columns_path)
    n_train, n_val = s_train.shape[0], s_val.shape[0]

    def gen_example_data_dummy_df(yx, s, begin, cols):
        """
        Generate dummy dataframe

        Args:
            yx (ndarray): [Y,X]
            s (ndarray): vector for site information
            begin (int): Begin RID
            cols (list): Columns for generated dataframe

        Returns:
            _type_: _description_
        """
        nb_subs = yx.shape[0]
        rids = np.array([i for i in range(begin, begin + nb_subs)]).reshape(
            (-1, 1))
        dates = np.array(['2022-01-01' for _ in range(nb_subs)]).reshape((-1,
                                                                          1))
        rng = np.random.default_rng(seed=0)
        dxs = (rng.integers(low=0, high=3, size=nb_subs)).reshape((-1, 1))
        dxs = dxs.astype(int)
        rng = np.random.default_rng(seed=0)
        ages = rng.random((nb_subs, 1))
        rng = np.random.default_rng(seed=0)
        sexs = (rng.integers(low=0, high=2, size=nb_subs)).reshape((-1, 1))
        sexs = sexs.astype(int)
        s = s.astype(int)
        array = np.concatenate((rids, dates, dxs, s, ages, sexs, yx), axis=1)
        df = pd.DataFrame(data=array, columns=cols)
        return df

    # for train
    train_df = gen_example_data_dummy_df(train_yx, s_train, 0, cols)
    # for val
    val_df = gen_example_data_dummy_df(val_yx, s_val, n_train, cols)
    # for test
    test_df = gen_example_data_dummy_df(test_yx, s_test, n_train + n_val, cols)

    return train_df, val_df, test_df


def gen_example_data(args):
    """
    Generate toy example data

    Args:
        args (tuple): Parameters
    """
    data_path = os.path.join(args.working_dir, 'data')
    raw_data_path = os.path.join(data_path, 'splits', args.dataset_pair)
    harm_input_path = os.path.join(data_path, 'harm_input', args.dataset_pair)
    goalDNN_input_path = os.path.join(data_path, 'goalDNN_input',
                                      args.dataset_pair)
    harm_output_path = os.path.join(data_path, 'harm_output',
                                    args.dataset_pair)
    create_folder(os.path.join(raw_data_path, '0'))
    create_folder(harm_input_path)
    create_folder(harm_output_path)
    create_folder(goalDNN_input_path)
    # step1, use sklearn to generate two heterogeneous datatases
    x_site0, y_site0, x_site1, y_site1 = gen_example_data_step1(args)
    # step2, split train/val/test
    train_yx, s_train, val_yx, s_val, test_yx, s_test =\
        gen_example_data_step2(args, x_site0, y_site0, x_site1, y_site1)
    # step3, generate dataframe as format of real ADNI data
    train_df, val_df, test_df = gen_example_data_step3(
        train_yx, s_train, val_yx, s_val, test_yx, s_test)
    # step4, save as csv format
    save_df(train_df,
            os.path.join(raw_data_path, '0', 'unmatch2match_train.csv'))
    save_df(val_df, os.path.join(raw_data_path, '0', 'unmatch2match_val.csv'))
    save_df(test_df, os.path.join(raw_data_path, '0',
                                  'unmatch2match_test.csv'))
    # step5, call utils.input_generation to generate data
    # 1. goalDNN unharmonized input
    goalDNN_unharmed_input_gen(
        raw_data_path, goalDNN_input_path, args.dataset_pair, nb_folds=1)
    # 2. VAE input
    VAE_input_gen(raw_data_path, harm_input_path, nb_folds=1)


def example_wrapper(args):
    """
    Wrapper function for running example/unit test of An2022_gcVAE project

    Args:
        args (tuple): Parameters
    """
    models = ['goalDNN', 'cVAE', 'gcVAE', 'XGBoost']
    stages = ['prepare', 'train', 'predict', 'print']

    assert args.model in models, 'Wrong input arguments for args.model'
    assert args.stage in stages, 'Wrong input arguments for args.stage'

    if args.unittest:
        args.working_dir = os.path.join(global_config.root_path, 'unit_tests')

    if args.stage == 'prepare':
        # generate toy example data
        gen_example_data(args)
    elif args.stage == 'train':
        if args.model == 'goalDNN':
            train_goaldnn_args = train_goalDNN_args_parser()
            train_goaldnn_args.data_path =\
                os.path.join(args.working_dir, 'data', 'goalDNN_input',
                             args.dataset_pair, '0')
            train_goaldnn_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints', 'eval_model',
                             'goalDNN', args.dataset_pair, '0')
            train_goaldnn_args.isSaving = True
            train_goaldnn_args.cpu = True
            train_goaldnn_args.epochs = 20
            train_goaldnn_args.lr = 1e-3
            train_goaldnn_args.lambda_dx = 1e-9
            train_goaldnn_args.lambda_mmse = 1
            train_goaldnn(train_goaldnn_args)
        elif args.model == 'cVAE':
            train_cvae_args = train_cvae_args_parser()
            train_cvae_args.data_path =\
                os.path.join(args.working_dir, 'data', 'harm_input',
                             args.dataset_pair, '0')
            train_cvae_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints', 'harm_model',
                             'cVAE', args.dataset_pair, '0')
            train_cvae_args.isSaving = True
            train_cvae_args.cpu = True
            train_cvae_args.epochs = 20
            train_cvae_args.lr = 1e-3
            train_cvae(train_cvae_args)
        elif args.model == 'gcVAE':
            train_gcvae_args = train_gcVAE_args_parser()
            train_gcvae_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints', 'harm_model',
                             'gcVAE', args.dataset_pair, '0')
            train_gcvae_args.goalDNN_input_path =\
                os.path.join(args.working_dir, 'data', 'goalDNN_input',
                             args.dataset_pair, '0')
            train_gcvae_args.harm_input_path =\
                os.path.join(args.working_dir, 'data', 'harm_input',
                             args.dataset_pair, '0')
            train_gcvae_args.cVAE_model =\
                os.path.join(args.working_dir, 'checkpoints', 'harm_model',
                             'cVAE', args.dataset_pair, '0', 'cVAE.pt')
            train_gcvae_args.goalDNN_model =\
                os.path.join(args.working_dir, 'checkpoints', 'eval_model',
                             'goalDNN', args.dataset_pair, '0', 'goalDNN.pt')
            train_gcvae_args.isSaving = True
            train_gcvae_args.cpu = True
            train_gcvae_args.epochs = 50
            train_gcvae_args.lr = 1e-4
            train_gcvae_args.lambda_dx = 1e-8
            train_gcvae_args.lambda_mmse = 1
            train_gcvae(train_gcvae_args)
        elif args.model == 'XGBoost':
            train_xgboost_args = train_XGBoost_args_parser()
            train_xgboost_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints',
                             'eval_model', 'XGBoost', args.dataset_pair)
            train_xgboost_args.output_path =\
                os.path.join(args.working_dir,
                             'results', 'XGBoost', args.dataset_pair)
            train_xgboost_args.nb_folds = 1
            train_xgboost_args.num_boost_rounds = 20
            train_xgboost_args.norm = False
            # for unharmonized
            train_xgboost_args.data_path =\
                os.path.join(args.working_dir, 'data',
                             'splits', args.dataset_pair)
            train_xgboost_args.train_name = 'unmatch2match_train'
            train_xgboost_args.val_name = 'unmatch2match_val'
            train_xgboost_args.test_name = 'unmatch2match_test'
            train_xgboost_args.save_suffix = 'unharm'
            train_XGBoost(train_xgboost_args)
            # for cVAE harmonization
            train_xgboost_args.data_path =\
                os.path.join(args.working_dir, 'data', 'harm_output',
                             args.dataset_pair, 'cVAE')
            train_xgboost_args.train_name = 'unmatch2match_train-intermediate'
            train_xgboost_args.val_name = 'unmatch2match_val-intermediate'
            train_xgboost_args.test_name = 'unmatch2match_test-intermediate'
            train_xgboost_args.save_suffix = 'cVAE'
            train_XGBoost(train_xgboost_args)
            # for gcVAE harmonization
            train_xgboost_args.data_path =\
                os.path.join(args.working_dir, 'data', 'harm_output',
                             args.dataset_pair, 'gcVAE')
            train_xgboost_args.train_name = 'unmatch2match_train-intermediate'
            train_xgboost_args.val_name = 'unmatch2match_val-intermediate'
            train_xgboost_args.test_name = 'unmatch2match_test-intermediate'
            train_xgboost_args.save_suffix = 'gcVAE'
            train_XGBoost(train_xgboost_args)

        else:
            raise myException('Unexpected --model')
    elif args.stage == 'predict':
        # making prediction using trained models
        if args.model == 'cVAE':
            predict_cvae_args = predict_cvae_args_parser()
            predict_cvae_args.raw_data_path =\
                os.path.join(
                    args.working_dir, 'data', 'splits')
            predict_cvae_args.harm_input_path =\
                os.path.join(args.working_dir, 'data', 'harm_input')
            predict_cvae_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints')
            predict_cvae_args.harm_output_path =\
                os.path.join(args.working_dir, 'data', 'harm_output')
            predict_cvae_args.nb_folds = 1
            predict_cvae(predict_cvae_args)
            # generate cVAE harmonized input for goalDNN
            test_name = 'unmatch2match_test-map2ADNI'
            goalDNN_harmed_input_gen(
                os.path.join(args.working_dir, 'data', 'harm_output',
                             args.dataset_pair, 'cVAE'),
                os.path.join(args.working_dir, 'data', 'goalDNN_input',
                             args.dataset_pair),
                args.dataset_pair,
                test_name,
                'cVAE',
                nb_folds=1)
        elif args.model == 'gcVAE':
            predict_gcvae_args = predict_gcvae_args_parser()
            predict_gcvae_args.raw_data_path =\
                os.path.join(args.working_dir, 'data', 'splits')
            predict_gcvae_args.harm_input_path =\
                os.path.join(args.working_dir, 'data', 'harm_input')
            predict_gcvae_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints')
            predict_gcvae_args.harm_output_path =\
                os.path.join(args.working_dir, 'data', 'harm_output')
            predict_gcvae_args.nb_folds = 1
            predict_gcvae(predict_gcvae_args)
            # generate gcVAE harmonized input for goalDNN
            test_name = 'unmatch2match_test-map2ADNI'
            goalDNN_harmed_input_gen(
                os.path.join(args.working_dir, 'data', 'harm_output',
                             args.dataset_pair, 'gcVAE'),
                os.path.join(args.working_dir, 'data', 'goalDNN_input',
                             args.dataset_pair),
                args.dataset_pair,
                test_name,
                'gcVAE',
                nb_folds=1)
        elif args.model == 'goalDNN':
            create_folder(
                os.path.join(args.working_dir, 'results', 'goalDNN',
                             args.dataset_pair))
            # making downstream application predction using trained goalDNN
            predict_goaldnn_args = predict_goalDNN_args_parser()
            predict_goaldnn_args.data_path =\
                os.path.join(args.working_dir, 'data', 'goalDNN_input',
                             args.dataset_pair)
            predict_goaldnn_args.checkpoint_path =\
                os.path.join(args.working_dir, 'checkpoints', 'eval_model',
                             'goalDNN', args.dataset_pair)
            predict_goaldnn_args.save_path =\
                os.path.join(args.working_dir,
                             'results', 'goalDNN', args.dataset_pair)
            predict_goaldnn_args.exp = 'unmatch2match'
            predict_goaldnn_args.nb_folds = 1
            predict_goaldnn_args.testsets = [
                'AIBL_unharm', 'AIBLcVAE_harm', 'AIBLgcVAE_harm'
            ]
            predict_goaldnn(predict_goaldnn_args)
            # evaluation
            eval_args = eval_goalDNN_args_parser()
            eval_args.save_path = os.path.join(args.working_dir, 'results',
                                               'goalDNN', args.dataset_pair)
            eval_args.nb_folds = 1
            eval_args.testsets = [
                'AIBL_unharm', 'AIBLcVAE_harm', 'AIBLgcVAE_harm'
            ]
            evaluation(eval_args.save_path, eval_args.testsets,
                       eval_args.nb_folds)
            for f in eval_args.testsets:
                f_name = 'DX_pred_result_' + f + '.txt'
                f_path = os.path.join(eval_args.save_path, f_name)
                os.remove(f_path)

        else:
            raise myException('Unexpected --model')
    else:
        # print results
        results_path = os.path.join(args.working_dir, 'results')
        mmse_models = ['AIBL_unharm', 'AIBLcVAE_harm', 'AIBLgcVAE_harm']
        mmse_prefix = 'MMSE_pred_result_'
        mmse_loger = []
        for mmmse_model in mmse_models:
            f_path = os.path.join(results_path, 'goalDNN', args.dataset_pair,
                                  mmse_prefix + mmmse_model + '.txt')
            mmse_result = float(((txt2list(f_path)[-1]).split('_'))[0])
            mmse_result = round(mmse_result, 4)
            mmse_loger.append(mmse_result)
        site_models = ['unharm', 'cVAE', 'gcVAE']
        site_prefix = 'site_pred_'
        site_loger = []
        for site_model in site_models:
            f_path = os.path.join(results_path, 'XGBoost', args.dataset_pair,
                                  site_prefix + site_model + '.txt')
            site_result = float(((txt2list(f_path)[-1]).split('_'))[0])
            site_result = round(site_result, 4)
            site_loger.append(site_result)
        # print
        print('>>>>>>> Without Harmonization <<<<<<<<')
        print('MMSE prediction MAE(lower=better):', mmse_loger[0])
        print('Dataset prediction accuracy(lower=better):', site_loger[0])
        print('>>>>>>> cVAE Harmonization <<<<<<<<')
        print('MMSE prediction MAE(lower=better):', mmse_loger[1])
        print('Dataset prediction accuracy(lower=better):', site_loger[1])
        print('>>>>>>> gcVAE Harmonization <<<<<<<<')
        print('MMSE prediction MAE(lower=better):', mmse_loger[2])
        print('Dataset prediction accuracy(lower=better):', site_loger[2])


if __name__ == '__main__':
    example_wrapper(example_args_parser())
