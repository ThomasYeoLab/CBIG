#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import random
import argparse
import torch
import numpy as np
from torch.nn import L1Loss, CrossEntropyLoss
from config import global_config
from utils.nn_misc import train_dataloader
from utils.metrics import subject_acc, subject_mae
from utils.misc import load_pkl, txt2list, list2txt, \
    one_hot, create_folder

loss_mae = L1Loss(reduction='mean')
loss_ce = CrossEntropyLoss(reduction='mean')


def train_gcVAE_args_parser():
    """
    Parameters for training gcVAE model
    """
    parser = argparse.ArgumentParser(prog='TraingcVAEArgs')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--GPU', type=int, default=-1)
    parser.add_argument('--node', '-n', type=int, default=0)
    parser.add_argument('--set', '-s', type=int, default=0)
    parser.add_argument('--isSaving', action='store_true', default=False)
    parser.add_argument('--cpu', action='store_true', default=False)
    parser.add_argument('--exp', type=str, default="unmatch2match")
    parser.add_argument('--dataset_pair', '-p', type=str, default='ADNI-AIBL')
    parser.add_argument('--cVAE_model', type=str, default='/')
    parser.add_argument('--goalDNN_model', type=str, default='/')
    parser.add_argument('--harm_input_path', type=str, default='/')
    parser.add_argument('--goalDNN_input_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    # training parameters
    parser.add_argument('--in_dim', type=int, default=108)
    parser.add_argument('--nb_sites', type=int, default=2)
    parser.add_argument('--epochs', type=int, default=1000)
    parser.add_argument('--step', type=int, default=0)
    parser.add_argument('--weight_decay', type=float, default=1e-6)
    parser.add_argument('--batch_size', type=int, default=128)
    # hyperparameters
    parser.add_argument('--lr', type=float, default=1e-6)
    parser.add_argument('--lr_step', type=int, default=100)
    parser.add_argument('--lambda_dx', type=float, default=1.0)
    parser.add_argument('--lambda_mmse', type=float, default=1.0)

    args, _ = parser.parse_known_args()

    return args


def train_1epoch(args, cVAE, goalDNN, dataloader, optimizer, task_mean_std,
                 harm_mean_std, mmse_std, device):
    """
    Train gcVAE model for one epoch

    Args:
        args (tuple): Parameters
        cVAE (class VAE): cVAE model
        goalDNN (class goalDNN): goalDNN model
        dataloader (class Dataloader): Training dataloader
        optimizer (class Adam): Adam optimizer
        task_mean_std (tuple): Mean and std from training goalDNN
        harm_mean_std (tuple): Mean and std from training cVAE
        mmse_std (tensor): Std of MMSE
        device (class Device): Device to train model
    """
    for i, batch_data in enumerate(dataloader):
        seed = args.step * len(dataloader) + i
        batch_ROIs = batch_data[:, :args.in_dim]
        batch_SITEs = batch_data[:, args.in_dim:args.in_dim + args.nb_sites]
        batch_MMSEs = \
            batch_data[:,
                       args.in_dim+args.nb_sites:args.in_dim+args.nb_sites+1]
        batch_DXs = batch_data[:, args.in_dim + args.nb_sites + 1:]
        # cVAE makes prediction
        cVAE.zero_grad()
        batch_ROIs = batch_ROIs.to(device)
        batch_SITEs = batch_SITEs.to(device)
        batch_MMSEs = batch_MMSEs.to(device)
        batch_DXs = batch_DXs.long()
        batch_DXs = torch.reshape(batch_DXs, (batch_DXs.shape[0], ))
        batch_DXs = batch_DXs.to(device)
        [ROIs_hat, _, _, _] = cVAE(batch_ROIs, batch_SITEs, seed)
        # denoramlization
        ROIs_hat = (ROIs_hat * harm_mean_std[1]) + harm_mean_std[0]
        ROIs_hat[ROIs_hat <= 0] = 0
        # renormalization
        ROIs_hat = (ROIs_hat - task_mean_std[0]) / task_mean_std[1]
        # goalDNN makes prediction
        goalDNN.zero_grad()
        [MMSE_pred, DX_pred] = goalDNN(ROIs_hat)
        # calculate loss
        mae_loss = loss_mae(MMSE_pred, batch_MMSEs) * mmse_std
        crossentropy_loss = loss_ce(DX_pred, batch_DXs)
        loss = args.lambda_mmse * mae_loss + args.lambda_dx * crossentropy_loss
        # backprogagation and update parameters
        loss.backward()
        optimizer.step()

    return cVAE, optimizer


def train(args):
    """
    Train gcVAE model using given hyper-parameters

    Args:
        args (tuple): Parameters
    """
    # set random seed
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    torch.cuda.manual_seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True
    # Set device to train on
    # By default, we use GPU to train model
    if args.GPU >= 0:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.GPU)
    else:
        # let gpuQ to select which GPU to train
        pass
    if args.cpu:
        device = torch.device('cpu')
    else:
        device = torch.device('cuda')
    create_folder(args.checkpoint_path, isOverwrite=False)
    # load model
    cVAE = torch.load(args.cVAE_model)
    goalDNN = torch.load(args.goalDNN_model)
    cVAE = cVAE.to(device)
    goalDNN = goalDNN.to(device)
    # load data
    ROIs = txt2list(global_config.ROI_features_path)
    train_pkl = load_pkl(os.path.join(args.harm_input_path, 'train.pkl'))
    val_pkl = load_pkl(os.path.join(args.harm_input_path, 'val_gcVAE.pkl'))
    # load harm mean and std
    harm_mean = train_pkl['mean'][ROIs].values
    harm_std = train_pkl['std'][ROIs].values
    harm_mmse_mean = train_pkl['mean'][['MMSE']].values
    harm_mmse_std = train_pkl['std'][['MMSE']].values
    # load task mean and std
    task_pkl = load_pkl(os.path.join(args.goalDNN_input_path, 'train.pkl'))
    task_mean = task_pkl['mean'][ROIs].values
    task_std = task_pkl['std'][ROIs].values
    if args.exp == 'unmatch2match':
        task_mmse_mean = train_pkl['mean'][['MMSE']].values
        task_mmse_std = train_pkl['std'][['MMSE']].values
    else:
        task_mmse_mean = task_pkl['mean'][['MMSE']].values
        task_mmse_std = task_pkl['std'][['MMSE']].values

    # training ROIs, MMSEs and DXs
    train_ROIs = train_pkl['data'][ROIs].values
    train_site_onehot = one_hot(train_pkl['data']['SITE'].values)
    train_MMSEs = train_pkl['data'][['MMSE']].values
    train_MMSEs = \
        (train_MMSEs * harm_mmse_std +
         harm_mmse_mean - task_mmse_mean) / task_mmse_std
    train_DXs = train_pkl['data'][['DX']].values
    train = np.concatenate(
        (train_ROIs, train_site_onehot, train_MMSEs, train_DXs), axis=1)
    # we only select nonADNI data for finetune !!#
    train = train[train[:, args.in_dim] == 0]
    # map nonADNI data to ADNI site
    train[:, args.in_dim] = 1
    train[:, args.in_dim + 1] = 0
    train = torch.tensor(train).float()

    # valdation ROIs, MMSEs and DXs
    val_ROIs = val_pkl['data'][ROIs].values
    val_index = np.where(val_pkl['data']['SITE'].values == 1)
    val_RIDs = val_pkl['RID'].values
    val_site_onehot = one_hot(val_pkl['data']['SITE'].values)
    # get MMSE and DX
    val_MMSEs = val_pkl['data'][['MMSE']].values
    val_MMSEs = \
        (val_MMSEs * harm_mmse_std +
         harm_mmse_mean - task_mmse_mean) / task_mmse_std
    val_DXs = val_pkl['data'][['DX']].values
    # concatenate Site, MMSE, DX
    val = np.concatenate((val_ROIs, val_site_onehot, val_MMSEs, val_DXs),
                         axis=1)
    # we only select nonADNI data for finetune !!#
    val = val[val[:, args.in_dim] == 0]
    # map nonADNI data to ADNI site
    val[:, args.in_dim] = 1
    val[:, args.in_dim + 1] = 0
    val = torch.tensor(val).float()
    val = val.to(device)

    harm_mean = torch.tensor(harm_mean).float()
    harm_mean = harm_mean.to(device)
    harm_std = torch.tensor(harm_std).float()
    harm_std = harm_std.to(device)
    harm_mean_std = (harm_mean, harm_std)

    task_mean = torch.tensor(task_mean).float()
    task_mean = task_mean.to(device)
    task_std = torch.tensor(task_std).float()
    task_std = task_std.to(device)
    task_mean_std = (task_mean, task_std)

    mmse_std = torch.tensor(task_mmse_std).float()
    mmse_std = mmse_std.to(device)

    # dataloader
    dataloader = train_dataloader(train, args.batch_size)
    # set optimizer, only update cVAE model !!!
    # optimizer
    optimizer = torch.optim.Adam(
        params=cVAE.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
        betas=(0.9, 0.999),
        eps=1e-7,
        amsgrad=False)
    lr_scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[args.lr_step], gamma=0.1)
    # val logger
    best_val = 1e5
    best_model = None
    val_logger = []
    # begin training
    for epoch in range(args.epochs):
        args.step = epoch
        cVAE.train()
        goalDNN.eval()
        # begin train for 1 epoch
        cVAE, optimizer = train_1epoch(args, cVAE, goalDNN, dataloader,
                                       optimizer, task_mean_std, harm_mean_std,
                                       mmse_std, device)
        lr_scheduler.step()
        val_ROIs = val[:, :args.in_dim]
        val_SITEs = val[:, args.in_dim:args.in_dim + args.nb_sites]
        val_MMSEs = \
            val[:, args.in_dim+args.nb_sites:args.in_dim+args.nb_sites+1]
        val_DXs = val[:, args.in_dim + args.nb_sites + 1:]
        valROIs_hat_mean = torch.zeros_like(val_ROIs)
        for i in range(100):
            [val_ROIs_hat, _, _, _] = cVAE(val_ROIs, val_SITEs, i)
            # data postprocessing
            val_ROIs_hat = \
                (val_ROIs_hat * harm_mean_std[1]) + harm_mean_std[0]
            val_ROIs_hat[val_ROIs_hat <= 0] = 0
            # renormalization
            val_ROIs_hat = \
                (val_ROIs_hat - task_mean_std[0]) / task_mean_std[1]
            valROIs_hat_mean += val_ROIs_hat / 100
        # valROIs_hat_mean = torch.nan_to_num(valROIs_hat_mean) # incase
        # goalDNN makes prediction
        [val_MMSEs_pred, val_DXs_pred] = goalDNN(valROIs_hat_mean)
        _, valMMSEMAE = subject_mae(val_MMSEs_pred, val_MMSEs,
                                    val_RIDs[val_index])
        valMMSEMAE *= mmse_std
        _, valDXAcc = subject_acc(val_DXs_pred, val_DXs, val_RIDs[val_index])
        if valMMSEMAE / 2 - valDXAcc < best_val:
            best_val = valMMSEMAE / 2 - valDXAcc
            best_model = cVAE
            log = \
                str(epoch) + '|' + str(best_val.data.cpu().numpy()) + \
                '|' + str(valMMSEMAE.data.cpu().numpy()) + \
                '|' + str(valDXAcc)
            val_logger.append(log)
    # end training
    if args.isSaving:
        torch.save(best_model, os.path.join(args.checkpoint_path, 'gcVAE2.pt'))
        val_logger_savename = str(
            'val_' + str(args.lr) + '_' + str(args.lr_step) + '_' +
            str(args.lambda_dx) + '_' + str(args.lambda_mmse) + '2.txt')
        list2txt(val_logger,
                 os.path.join(args.checkpoint_path, val_logger_savename))
    else:
        val_logger_savename = str(
            'val_' + str(args.lr) + '_' + str(args.lr_step) + '_' +
            str(args.lambda_dx) + '_' + str(args.lambda_mmse) + '.txt')
        list2txt(val_logger,
                 os.path.join(args.checkpoint_path, val_logger_savename))


if __name__ == '__main__':
    train(train_gcVAE_args_parser())
