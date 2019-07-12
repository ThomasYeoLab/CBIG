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
from pathlib import PurePath

import keras
from keras import backend as K
from keras.optimizers import SGD, Adam

from config import config
from CBIG_model_keras import gcnn
from CBIG_mics import mics_log, mics_graph_matrix, mics_eval


def train(args):
    '''main function for GCNN network ans sex

    Args:
        args: args from command line

    Returns:
        None
    '''

    t_overall = time.time()
    print('\nCBIG GCNN for UK Biobank and sex with argument: ' + str(args))

    # set seed (Keras seed does not work well to generate same prediction)
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

    # load data
    npz = PurePath(args.path_data, 'data_gcnn_sex.npz').as_posix()
    npz = np.load(npz)
    input_x = npz['input_x']
    input_y_original = np.squeeze(npz['input_y'])
    test_mask = npz['test_mask']
    valid_mask = npz['valid_mask']
    train_mask = npz['train_mask']

    # Convert y for Keras
    input_y = keras.utils.to_categorical(input_y_original, 2)

    # Define graph parameters
    runs = args.runs  # numbers of ensemble runs
    epochs = args.epochs  # numbers of epochs per run
    num_subject = args.num_subject
    test_index = np.nonzero(test_mask)[0]  # index to get test result
    graph_matrix, support = mics_graph_matrix(
        num_subject, args.graph_folder, args.graph_setup,
        args.graph_estimation, args.chebyshev_degree)

    # initialization of result record
    tra_los_record = np.zeros((runs, epochs))
    tra_auc_record = np.zeros((runs, epochs))
    val_auc_record = np.zeros((runs, epochs))
    tes_auc_record = np.zeros((runs, epochs))
    tes_res_record = np.zeros((runs, epochs, np.sum(test_mask), 2))
    final_original = None

    # Code running - with multiple ensemble runs
    for run in range(runs):
        K.clear_session()
        print('Run:', run + 1)

        model = gcnn(
            input_x.shape[1],
            args.dropout,
            args.n_l1,
            graph_matrix,
            support,
            args.l2_regularizer,
            for_sex=True)

        if args.optimizer == 'SGD':
            optimizer = SGD(
                lr=args.lr,
                momentum=args.momentum,
                decay=args.lr_decay,
                nesterov=False)
        elif args.optimizer == 'Adam':
            optimizer = Adam(lr=args.lr)
        else:
            raise Exception('Invalid optimizer type.')
        model.compile(
            loss='categorical_crossentropy',
            optimizer=optimizer,
            metrics=['accuracy'])

        # run each epoch
        for epoch in range(epochs):

            h = model.fit(
                input_x,
                input_y,
                sample_weight=train_mask,
                batch_size=num_subject,
                epochs=1,
                shuffle=False,
                verbose=0)
            preds = model.predict(input_x, batch_size=num_subject)

            val_auc, tes_auc, tra_auc = mics_eval(preds, input_y, train_mask,
                                                  valid_mask, test_mask)

            # log result
            tra_los_record[run, epoch] = h.history['loss'][0]
            tra_auc_record[run, epoch] = tra_auc
            val_auc_record[run, epoch] = val_auc
            tes_auc_record[run, epoch] = tes_auc
            tes_real = np.squeeze(input_y[test_index, :])
            if final_original is not None:
                assert np.array_equal(final_original, tes_real)
            else:
                final_original = tes_real
            tes_res_record[run, epoch, :, :] = np.squeeze(preds[test_index, :])

    log_args = {
        'tra_los_record': tra_los_record,
        # 'val_los_record': val_los_record,
        # 'tes_los_record': tes_los_record,
        'tra_auc_record': tra_auc_record,
        'val_auc_record': val_auc_record,
        'tes_auc_record': tes_auc_record,
        'tes_res_record': tes_res_record,
        'final_original': final_original
    }

    mics_log('gcnn', args.out_path, index=args.index, **log_args)
    print("time spent: {:.4f}".format(time.time() - t_overall))

    return


def get_args():
    '''function to get args from command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--path_data', type=str, default=config.UKBB_INTER_DIR)
    parser.add_argument('--out_path', '-o', type=str, default=config.OUT_PATH)
    parser.add_argument(
        '--graph_folder', type=str, default=config.GRAPH_FOLDER)
    parser.add_argument('--seed', type=int, default=config.RAMDOM_SEED)
    parser.add_argument('--runs', type=int, default=config.UKBB_RUNS)
    parser.add_argument(
        '--num_subject', type=int, default=config.UKBB_NUM_SUBJECT)
    parser.add_argument('--epochs', type=int, default=config.UKBB_EPOCHS_GCNN)
    parser.add_argument('--gpu', type=int, default=0)

    # hyperparameter
    parser.add_argument('--index', type=int, default=None)
    parser.add_argument(
        '--graph_setup', default='8868Subject_corr_option_3_param_5')
    parser.add_argument('--graph_estimation', default='chebyshev')
    parser.add_argument('--chebyshev_degree', type=int, default=1)
    parser.add_argument('--optimizer', default='Adam')
    parser.add_argument('--lr', type=float, default=0.0003)
    parser.add_argument('--momentum', type=float, default=0.9)
    parser.add_argument('--lr_decay', type=float, default=1e-6)
    parser.add_argument('--l2_regularizer', type=float, default=2e-5)
    parser.add_argument('--dropout', type=float, default=0.5)
    parser.add_argument('--n_l1', type=int, default=6)

    return parser.parse_args()


if __name__ == '__main__':
    train(get_args())
