#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import torch
import argparse
import numpy as np
from torch.nn import L1Loss, CrossEntropyLoss
from config import global_config
from model.goalDNN import goalDNN
from utils.nn_misc import train_dataloader, extract_goaldnn_input
from utils.metrics import subject_acc, subject_mae
from utils.misc import txt2list, load_pkl, create_folder

# loss functions used for training goalDNN model
loss_mae = L1Loss(reduction='mean')
loss_ce = CrossEntropyLoss(reduction='mean')


def train_goalDNN_args_parser():
    """
    Parameters for training goalDNN model
    """
    parser = argparse.ArgumentParser(prog='TraingoalDNNArgs')
    # general parameters
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--GPU', type=int, default=-1)
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--isSaving', action='store_true', default=False)
    parser.add_argument('--cpu', action='store_true', default=False)
    # training parameters
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--step', type=int, default=0)
    parser.add_argument('--batch_size', type=int, default=128)
    parser.add_argument('--in_dim', type=int, default=108)
    parser.add_argument('--nb_category', type=int, default=3)
    parser.add_argument('--nb_measures', type=int, default=1)
    parser.add_argument(
        '--hidden_dims', type=list, default=[512, 256, 128, 64, 32])

    parser.add_argument('--drop_out', type=float, default=0.1)
    parser.add_argument('--weight_decay', type=float, default=1e-4)
    # hyperparameters
    parser.add_argument('--lr', type=float, default=0.000174138)
    parser.add_argument('--lambda_dx', type=float, default=0.9)
    parser.add_argument('--lambda_mmse', type=float, default=0.4)
    parser.add_argument('--lr_step', type=int, default=95)
    parser.add_argument('--h1', type=int, default=512)
    parser.add_argument('--h2', type=int, default=64)
    parser.add_argument('--h3', type=int, default=512)
    parser.add_argument('--h4', type=int, default=512)
    parser.add_argument('--h5', type=int, default=512)
    parser.add_argument('--nb_layers', type=int, default=5)

    args, _ = parser.parse_known_args()
    return args


def train_1epoch(args, dataloader, model, optimizer, device, mmse_std):
    """
    Train goalDNN model for one epoch

    Args:
        args (tuple): Parameters
        dataloader (class Dataloader): Training dataloader
        model (class VAE): cVAE model
        optimizer (class Adam): Adam optimizer
        device (class device): Device the model training on
        mmse_std (tensor): std of MMSE
    """
    for _, batch_data in enumerate(dataloader):
        batchROIs = batch_data[:, :args.in_dim].float()
        batchMMSEs = \
            batch_data[:, args.in_dim:args.in_dim+args.nb_measures].float()
        batchDXs = batch_data[:, args.in_dim + args.nb_measures:].long()
        batchDXs = torch.reshape(batchDXs, (batchDXs.shape[0], ))
        # move to GPU
        batchROIs = batchROIs.to(device)
        batchMMSEs = batchMMSEs.to(device)
        batchDXs = batchDXs.to(device)
        model.zero_grad()
        # forward pass
        [MMSE_pred, DX_pred] = model(batchROIs)
        # calculate losses
        mae_loss = loss_mae(MMSE_pred, batchMMSEs) * mmse_std
        crossentropy_loss = \
            loss_ce(DX_pred, batchDXs)
        loss = args.lambda_mmse * mae_loss + args.lambda_dx * crossentropy_loss
        # backward & update parameters
        loss.backward()
        optimizer.step()
    return model, optimizer


def train(args):
    """
    Wrapper function for training goalDNN model

    Args:
        args (tuple): Parameters
    """
    if args.isSaving:
        create_folder(args.checkpoint_path)
    hid_dims = []
    hid_dims.append(args.h1)
    hid_dims.append(args.h2)
    hid_dims.append(args.h3)
    hid_dims.append(args.h4)
    hid_dims.append(args.h5)
    args.hidden_dims = hid_dims[:args.nb_layers]
    # Set random seed
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    # Set device to train on
    # By default, we use GPU to train model
    if args.GPU >= 0:
        pass
    else:
        # let gpuQ to select which GPU to train
        pass
    if args.cpu:
        device = torch.device('cpu')
    else:
        device = torch.device('cuda')
    # load data
    # read features of training data ==> 108 brain ROI volumes
    ROI_features = txt2list(global_config.ROI_features_path)
    # load training data
    train_pkl = load_pkl(os.path.join(args.data_path, 'train.pkl'))
    val_pkl = load_pkl(os.path.join(args.data_path, 'val.pkl'))
    # extract data for model training
    trainROIs, trainMMSEs, trainDXs, _ = \
        extract_goaldnn_input(train_pkl, ROI_features)
    valROIs, valMMSEs, valDXs, val_index = \
        extract_goaldnn_input(val_pkl, ROI_features)
    valROIs = valROIs.float()
    valROIs = valROIs.to(device)
    valMMSEs = valMMSEs.float()
    valMMSEs = valMMSEs.to(device)
    valDXs = valDXs.long()
    valDXs = valDXs.to(device)

    # build model
    model = goalDNN(
        in_dim=args.in_dim,
        nb_category=args.nb_category,
        nb_measures=args.nb_measures,
        p_dropout=args.drop_out,
        hidden_dims=args.hidden_dims)
    model = model.to(device)
    # create dataloader
    train_data = torch.cat((trainROIs, trainMMSEs, trainDXs), dim=1)
    dataloader = train_dataloader(train_data, args.batch_size)
    # optimizer
    optimizer = torch.optim.Adam(
        params=model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
        betas=(0.9, 0.999),
        eps=1e-7,
        amsgrad=False)
    lr_scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[args.lr_step], gamma=0.1)
    # best_val
    best_model = None
    best_val = 1e5
    best_valMMSEMAE = 1e5
    best_valDXAcc = 0
    # begin training
    for epoch in range(args.epochs):
        args.step = epoch
        # train 1 epoch
        model.train()
        model, optimizer = \
            train_1epoch(args, dataloader, model,
                         optimizer, device, train_pkl['std']['MMSE'])
        lr_scheduler.step()
        # evaluate model performance on validation set
        model.eval()
        [valMMSEs_pred, valDXs_pred] = model(valROIs)
        _, valMMSEMAE = subject_mae(valMMSEs_pred, valMMSEs,
                                    val_pkl['RID'].values[val_index])
        valMMSEMAE *= train_pkl['std']['MMSE']
        _, valDXAcc = subject_acc(valDXs_pred, valDXs,
                                  val_pkl['RID'].values[val_index])
        if valMMSEMAE / 2 - valDXAcc < best_val:
            best_val = valMMSEMAE / 2 - valDXAcc
            best_valMMSEMAE = valMMSEMAE
            best_valDXAcc = valDXAcc
            best_model = model

    if args.isSaving:
        torch.save(best_model, os.path.join(args.checkpoint_path,
                                            'goalDNN.pt'))
    else:
        # output via stdout
        print(str(best_valMMSEMAE) + ',' + str(best_valDXAcc))


if __name__ == '__main__':
    train(train_goalDNN_args_parser())
