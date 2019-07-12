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
from CBIG_model_keras import gcnn
from CBIG_mics import mics_hcp_log, mics_graph_matrix, mics_z_norm_gcnn
from CBIG_mics import mics_train_valid_mask_split, mics_hcp_gcnn_eval


def train(args):
    '''main function for GCNN network

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()
    print('\nCBIG GCNN for HCP with argument: ' + str(args))

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
    graph_matrix, support = mics_graph_matrix(
        args.num_subjects, args.graph_folder, args.graph_setup,
        args.graph_estimation, args.chebyshev_degree)

    # initialization of result record
    tra_los_log = np.zeros((tes_folds, val_folds, epochs))
    val_los_log = np.zeros((tes_folds, val_folds, epochs))
    val_cor_log = np.zeros((tes_folds, val_folds, epochs, n_measure))
    tra_cor_log = np.zeros((tes_folds, epochs, n_measure))
    tes_cor_log = np.zeros((tes_folds, epochs, n_measure))
    tes_mae_log = np.zeros((tes_folds, epochs, n_measure))
    final_orig = np.zeros((args.num_subjects, n_measure))
    final_pred = np.zeros((epochs, args.num_subjects, n_measure))
    test_count = np.zeros((tes_folds)).astype(int)
    t_sigma = np.zeros((tes_folds, n_measure))

    cnt = 0

    # cross validation start
    for tes_fold in range(tes_folds):

        # load data
        npz = os.path.join(args.path_data, 'gcnn',
                           'data_fold' + str(tes_fold + 1) + '.npz')
        npz = np.load(npz)
        input_x = npz['input_x']
        input_y_original = npz['input_y']
        test_mask = npz['test_mask']
        train_valid_mask = npz['train_valid_mask']

        # inner loop
        for val_fold in range(val_folds + 1):
            K.clear_session()

            # log
            if val_fold == val_folds:
                print('Test with train+valid:', tes_fold + 1)
                [train_mask, _] = mics_train_valid_mask_split(train_valid_mask)
            else:
                print('Test:', tes_fold + 1, 'Valid', val_fold + 1)
                [train_mask, valid_mask] = mics_train_valid_mask_split(
                    train_valid_mask, fold=val_fold)

            # z normalize data based on training set
            input_y, t_v_sigma = mics_z_norm_gcnn(input_y_original, train_mask)
            if val_fold == val_folds:
                test_count[tes_fold] = np.sum(test_mask)
                final_orig[cnt:(cnt + test_count[tes_fold]), :] =\
                    input_y[np.nonzero(test_mask)[0], :]
                t_sigma[tes_fold, :] = t_v_sigma

            model = gcnn(
                input_x.shape[1],
                args.dropout,
                args.n_l1,
                graph_matrix,
                support,
                args.l2_regularizer,
                n_measure=n_measure)
            optimizer = SGD(
                lr=args.lr,
                momentum=args.momentum,
                decay=args.lr_decay,
                nesterov=False)
            model.compile(loss='mean_squared_error', optimizer=optimizer)

            # run epochs
            for e in range(epochs):
                h = model.fit(
                    input_x,
                    input_y,
                    sample_weight=train_mask,
                    batch_size=args.num_subjects,
                    epochs=1,
                    shuffle=False,
                    verbose=0)
                preds = model.predict(input_x, batch_size=args.num_subjects)

                # calculate correlation
                if val_fold == val_folds:
                    tes_cor, _, tes_mae, tra_cor = mics_hcp_gcnn_eval(
                        preds, input_y, test_mask, t_v_sigma, train_mask)
                    tes_cor_log[tes_fold, e, :] = tes_cor
                    tra_cor_log[tes_fold, e, :] = tra_cor
                    tes_mae_log[tes_fold, e, :] = tes_mae
                    final_pred[e, cnt:(cnt + test_count[tes_fold]), :] =\
                        preds[np.nonzero(test_mask)[0], :]
                else:
                    val_cor, val_los, _ = mics_hcp_gcnn_eval(
                        preds, input_y, valid_mask)
                    tra_los_log[tes_fold, val_fold, e] = h.history['loss'][0]
                    val_los_log[tes_fold, val_fold, e] = val_los
                    val_cor_log[tes_fold, val_fold, e, :] = val_cor

        cnt = cnt + test_count[tes_fold]

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
    mics_hcp_log('gcnn', args.out_path, **log_args)
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
        '--batch_size', type=int, default=config.HCP_BATCH_SIZE)
    parser.add_argument(
        '--graph_folder', type=str, default=config.GRAPH_FOLDER)
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--gpu', type=int, default=0)

    # hyperparameter
    parser.add_argument(
        '--graph_setup', default='953Subject_corr_option_3_param_5')
    parser.add_argument('--graph_estimation', default='chebyshev')
    parser.add_argument('--chebyshev_degree', type=int, default=5)
    parser.add_argument('--optimizer', default='SGD')
    parser.add_argument('--lr', type=float, default=0.01)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--lr_decay', type=float, default=1e-5)
    parser.add_argument('--l2_regularizer', type=float, default=8e-4)
    parser.add_argument('--dropout', type=float, default=0.3)
    parser.add_argument('--n_l1', type=int, default=256)

    return parser.parse_args()


if __name__ == '__main__':
    train(get_args())
