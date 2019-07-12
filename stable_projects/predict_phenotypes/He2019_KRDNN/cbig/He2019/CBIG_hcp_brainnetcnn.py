#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import random
import argparse
import numpy as np
from keras import backend as K
from keras.optimizers import SGD

from config import config
from CBIG_model_keras import hcp_brainnetcnn
from CBIG_mics import mics_hcp_log, mics_train_valid_split
from CBIG_mics import mics_z_norm_test, mics_hcp_infer


def train(args):
    '''main function for BrainNetCNN network

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()
    print('\nCBIG BrainNetCNN for HCP with argument: ' + str(args))

    # set seed (Keras seed does not generate repeatable prediction)
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

    epochs = args.epochs  # numbers of epochs
    tes_folds = args.folds  # numbers of test folds
    val_folds = tes_folds - 1  # numbers of inner loop folds
    n_measure = args.n_measure

    # initialization of result record
    tra_los_log = np.zeros((tes_folds, val_folds, epochs))
    val_los_log = np.zeros((tes_folds, val_folds, epochs))
    val_cor_log = np.zeros((tes_folds, val_folds, epochs, n_measure))
    tra_cor_log = np.zeros((tes_folds, epochs, n_measure))
    tes_cor_log = np.zeros((tes_folds, epochs, n_measure))
    tes_mae_log = np.zeros((tes_folds, epochs, n_measure))
    final_orig = np.zeros((args.num_subjects, n_measure))
    final_pred = np.zeros((epochs, args.num_subjects, n_measure))
    test_count = np.zeros((tes_folds))
    t_sigma = np.zeros((tes_folds, n_measure))

    cnt = 0

    # cross validation start
    for tes_fold in range(tes_folds):

        # load data
        npz = os.path.join(args.path_data, 'brainnetcnn',
                           'data_fold' + str(tes_fold + 1) + '.npz')
        npz = np.load(npz)
        test_x = np.expand_dims(npz['test_x'], axis=-1)
        test_y = npz['test_y']
        train_valid_x = npz['train_valid_x']
        train_valid_y = npz['train_valid_y']
        test_y, t_v_sigma = mics_z_norm_test(train_valid_y, test_y)
        test_count[tes_fold] = test_y.shape[0]
        final_orig[cnt:(cnt + test_y.shape[0]), :] = test_y
        t_sigma[tes_fold, :] = t_v_sigma

        # inner loop
        for val_fold in range(val_folds + 1):
            K.clear_session()

            # log
            if val_fold == val_folds:
                print('Test with train+valid:', tes_fold + 1)
                train_x, train_y = mics_train_valid_split(
                    train_valid_x, train_valid_y, is_bnc=True)
            else:
                print('Test fold:', tes_fold, 'valid fold:', val_fold)
                train_x, valid_x, train_y, valid_y = mics_train_valid_split(
                    train_valid_x, train_valid_y, fold=val_fold, is_bnc=True)

            # initialize model
            model = hcp_brainnetcnn(train_x.shape[1], n_measure, args.e2e,
                                    args.e2n, args.n2g, args.dropout,
                                    args.leaky_alpha)
            optimizer = SGD(
                lr=args.lr,
                momentum=args.momentum,
                decay=args.lr_decay,
                nesterov=False)
            model.compile(loss='mean_squared_error', optimizer=optimizer)

            # run epochs
            for e in range(epochs):
                # fit model
                if val_fold == val_folds:
                    h = model.fit(
                        train_x,
                        train_y,
                        epochs=1,
                        batch_size=args.batch_size,
                        verbose=0)
                else:
                    h = model.fit(
                        train_x,
                        train_y,
                        epochs=1,
                        batch_size=args.batch_size,
                        verbose=0,
                        validation_data=(valid_x, valid_y))

                # calculate correlation
                if val_fold == val_folds:
                    tes_cor, tes_mae, y_pred, tra_cor = mics_hcp_infer(
                        model, test_x, test_y, t_v_sigma, train_x, train_y)
                    tes_cor_log[tes_fold, e, :] = tes_cor
                    tra_cor_log[tes_fold, e, :] = tra_cor
                    tes_mae_log[tes_fold, e, :] = tes_mae
                    final_pred[e, cnt:(cnt + test_y.shape[0]), :] = y_pred
                else:
                    val_cor, _, _ = mics_hcp_infer(model, valid_x, valid_y,
                                                   t_v_sigma)
                    tra_los_log[tes_fold, val_fold, e] = h.history['loss'][0]
                    val_los_log[tes_fold, val_fold,
                                e] = h.history['val_loss'][0]
                    val_cor_log[tes_fold, val_fold, e, :] = val_cor

        cnt = cnt + test_y.shape[0]

    log_args = {
        'tra_los_log': tra_los_log,
        'val_los_log': val_los_log,
        'tra_cor_log': tra_cor_log,
        'val_cor_log': val_cor_log,
        'tes_cor_log': tes_cor_log,
        'tes_mae_log': tes_mae_log,
        't_sigma': t_sigma,
        'final_orig': final_orig,
        'final_pred': final_pred,
        'test_count': test_count
    }
    mics_hcp_log('brainnetcnn', args.out_path, **log_args)
    print("time spent: {:.4f}".format(time.time() - t_overall))

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--path_data', type=str, default=config.HCP_INTER_DIR)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument(
        '--num_subjects', type=int, default=config.HCP_NUM_SUBJECT)
    parser.add_argument('--out_path', '-o', default=config.OUT_PATH)
    parser.add_argument('--folds', type=int, default=config.HCP_NUM_FOLD)
    parser.add_argument('--n_measure', type=int, default=config.HCP_N_MEASURE)
    parser.add_argument(
        '--batch_size', type=int, default=config.UKBB_BATCH_SIZE)
    parser.add_argument('--epochs', type=int, default=150)
    parser.add_argument('--gpu', type=int, default=0)

    # hyperparameter
    parser.add_argument('--lr', type=float, default=0.01)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--lr_decay', type=float, default=1e-4)
    parser.add_argument('--leaky_alpha', type=float, default=0.1)
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--e2e', type=int, default=16)
    parser.add_argument('--e2n', type=int, default=128)
    parser.add_argument('--n2g', type=int, default=26)

    return parser.parse_args()


if __name__ == '__main__':
    train(get_args())
