#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import time
import random
import torch
import optuna
import numpy as np
import torch.utils.data
import torch.optim as optim
from optuna.samplers import TPESampler
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from config import config
from CBIG_model_pytorch import msenanloss, multi_task_dataset, dnn
from CBIG_misc import misc_z_norm, misc_infer_metric, demean_normalize


def objective(trial):
    '''Optuna optimization function to search for the best hyperparameters for the DNN network.

    This function is called by an Optuna trial to suggest hyperparameters and evaluate them by training a DNN model.

    Args:
        trial (optuna.trial.Trial): An Optuna trial object that suggests hyperparameters.

    Returns:
        float: The validation COD of the model trained with the suggested hyperparameters.
    '''

    t_overall = time.time()

    dropout = trial.suggest_float("dropout", 0, 0.6, step=0.1)
    n_l1 = trial.suggest_categorical("n_l1", [128, 256, 512, 1024])
    n_l2 = trial.suggest_categorical("n_l2", [128, 256, 512, 1024])
    n_l3 = trial.suggest_categorical("n_l3", [128, 256, 512, 1024])
    n_l4 = trial.suggest_categorical("n_l4", [128, 256, 512, 1024])
    n_layer = trial.suggest_int("n_layer", 0, 3, step=1)

    # load dataset for PyTorch
    batch_size = trial.suggest_categorical("batch_size", [128, 256, 512])
    dset_train = multi_task_dataset(x_train, y_train)
    trainloader = DataLoader(
        dset_train, batch_size=batch_size, shuffle=True, num_workers=0)
    dset_valid = multi_task_dataset(x_valid, y_valid)
    validLoader = DataLoader(
        dset_valid, batch_size=batch_size, shuffle=True, num_workers=0)

    runs = 1  # numbers of ensemble runs
    epochs = config.EPOCHS  # numbers of epochs per run

    # initialization of result record
    tra_los_record = np.zeros((runs, epochs))
    val_los_record = np.zeros((runs, epochs))
    val_cor_record = np.zeros((runs, epochs, n_phe))
    val_cod_record = np.zeros((runs, epochs, n_phe))
    val_mae_record = np.zeros((runs, epochs, n_phe))
    val_res_record = np.zeros((runs, epochs, x_valid.shape[0], n_phe))

    # Code running - with multiple ensemble runs
    for run in range(runs):
        # initialization of network
        net = dnn(
            x_train.shape[1],
            n_layer,
            n_l1,
            n_l2,
            n_l3,
            n_l4,
            dropout,
            output_size=n_phe)
        # print(net)
        net.to(device)

        # other components of network
        criterion = msenanloss

        optimizer = optim.SGD(
            net.parameters(),
            lr=trial.suggest_float("lr", 1e-7, 1e-1, log=True),
            momentum=0.9,
            weight_decay=trial.suggest_float(
                "weight_decay", 1e-7, 1e-1, log=True))

        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='min', factor=0.5, patience=25, verbose=True)

        # start epoch training
        for epoch in range(epochs):
            # training
            train_loss = 0.0
            net.train(True)
            for (x, y) in trainloader:
                x, y = x.to(device), y.to(device)
                optimizer.zero_grad()
                outputs = net(x)
                mask = torch.isnan(y)
                y.masked_fill_(mask, 0)
                loss = criterion(outputs, y, mask=mask)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
            tra_los_record[run, epoch] = train_loss / len(trainloader)

            net.train(False)

            # validation
            corr, cod, mae, loss, _, pred = misc_infer_metric(
                validLoader,
                net,
                criterion,
                device,
                t_sigma,
                output_size=n_phe,
                need_value=True)
            val_cor_record[run, epoch, :] = corr
            val_cod_record[run, epoch, :] = cod
            val_mae_record[run, epoch, :] = mae
            val_los_record[run, epoch] = loss
            val_res_record[run, epoch, :, :] = np.squeeze(pred)

            scheduler.step(loss)
            trial.report(np.mean(cod), epoch)
            if trial.should_prune():
                raise optuna.exceptions.TrialPruned()

    print("time spent: {:.4f}".format(time.time() - t_overall))
    temp = np.mean(val_cod_record[0, :, :], axis=1)
    temp = np.convolve(temp, np.ones(3, dtype=int), 'valid') / 3
    index = np.nanargmax(temp)
    index = index + 1
    print('\nBest validation at index: ', index)
    best_cod = np.mean(val_cod_record[0, index, :])
    return best_cod


if __name__ == '__main__':
    # set gpu number
    os.environ["CUDA_VISIBLE_DEVICES"] = str(0)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # set all the seed
    seed = config.RAMDOM_SEED
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    # load data
    npz = os.path.join(config.IN_DIR, 'UKBB', 'UKBB_dnn_input.npz')
    npz = np.load(npz)
    x_UKBB = npz['x_raw']
    x_UKBB[np.isnan(x_UKBB)] = 0
    y_UKBB = npz['y_raw']

    # subject-wise normalization for input functional connectivity
    x_UKBB = demean_normalize(x_UKBB)
    print(x_UKBB.shape, y_UKBB.shape)

    # split train and validation
    split_tra, split_val = train_test_split(
        np.arange(x_UKBB.shape[0]), test_size=0.2, random_state=seed)

    x_train = x_UKBB[split_tra, :]
    x_valid = x_UKBB[split_val, :]
    y_train = y_UKBB[split_tra, :]
    y_valid = y_UKBB[split_val, :]

    # z norm based on y_train
    y_train, y_valid, _, t_sigma = misc_z_norm(y_train, y_valid)
    n_phe = y_train.shape[1]

    sampler = TPESampler(seed=seed)
    study = optuna.create_study(sampler=sampler, direction='maximize')
    study.optimize(objective, n_trials=200, gc_after_trial=True)
    trial = study.best_trial
    print('Accuracy: {}'.format(trial.value))
    print("Best hyperparameters: {}".format(trial.params))
