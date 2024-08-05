#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import os
import torch
import random
import argparse
import numpy as np
from torch.nn import functional as F
from config import global_config
from model.VAE import VAE
from model.VAE_modules import Discriminator
from utils.misc import \
    load_pkl, txt2list, one_hot, list2txt, create_folder
from utils.nn_misc import \
    KLD_loss, kl_conditional_and_marg, discriminator_loss, train_dataloader


def train_covae_args_parser():
    """
    Parameters for training cVAE model
    """
    parser = argparse.ArgumentParser(prog='TraincoVAEArgs')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--GPU', type=int, default=1)
    parser.add_argument('--data_path', type=str, default='/')
    parser.add_argument('--checkpoint_path', type=str, default='/')
    parser.add_argument('--model_name', type=str, default='coVAE')
    parser.add_argument('--isSaving', action='store_true', default=False)
    # training parameters
    parser.add_argument('--epochs', type=int, default=1000)
    parser.add_argument('--adv_epochs', type=int, default=1)
    parser.add_argument('--step', type=int, default=0)
    parser.add_argument('--batch_size', type=int, default=512)
    parser.add_argument('--nb_classes', type=int, default=9)
    parser.add_argument('--in_dim', type=int, default=108)
    parser.add_argument('--VAE_hidden_dims', type=list, default=[96, 64])
    parser.add_argument('--D_hidden_dims', type=list, default=[32, 32])
    parser.add_argument('--dis_lr', type=float, default=1e-2)
    parser.add_argument('--coff', type=float, default=0)
    # hyper-parameters
    parser.add_argument('--lr', type=float, default=1e-2)
    parser.add_argument('--drop_out', type=float, default=0.1)
    parser.add_argument('--alpha', type=float, default=0.001)
    parser.add_argument('--lambda_', type=float, default=0.001)
    parser.add_argument('--gamma', type=float, default=1.)
    parser.add_argument('--lr_step', type=int, default=100)
    parser.add_argument('--latent_dim', type=int, default=32)
    parser.add_argument('--h1', type=int, default=512)
    parser.add_argument('--h2', type=int, default=512)
    parser.add_argument('--h3', type=int, default=512)
    parser.add_argument('--h4', type=int, default=512)
    parser.add_argument('--nb_layers', type=int, default=3)

    train_args, _ = parser.parse_known_args()
    return train_args


def train_1epoch(args, dataloader, models, optimizers, device):
    """
    Train cVAE model for one epoch
    """
    net_G, net_D = models
    optimizer_G, optimizer_D = optimizers
    for i, batch_data in enumerate(dataloader):
        batch_x = batch_data[:, :args.in_dim].float()
        batch_y = batch_data[:, args.in_dim:].float()
        batch_x = batch_x.to(device)
        batch_y = batch_y.to(device)
        seed = args.step * len(dataloader) + i
        ####################################################
        # Train Discriminator k times (k=args.nb_adv_epochs)
        ####################################################
        for _ in range(args.adv_epochs):
            net_D.zero_grad()
            # train with real samples: batch_x
            real_labels = torch.ones((batch_x.shape[0]), device=device)
            real_pred = net_D(batch_x, args.coff)
            real_loss = discriminator_loss(real_pred, real_labels)
            real_loss.backward()
            # train with fake samples : reconstructed x_hat
            [x_hat, _, _, _] = net_G(batch_x, batch_y, seed)
            fake_labels = torch.zeros((x_hat.shape[0]), device=device)
            # using detach() to avoid accumulating gradient on Generator
            fake_pred = net_D(x_hat.detach(), args.coff)
            fake_loss = discriminator_loss(fake_pred, fake_labels)
            fake_loss.backward()
            optimizer_D.step()
        ####################################################
        # Train Generator (VAE)
        ####################################################
        net_G.zero_grad()
        [x_hat, x, z_mean, z_logvar] = net_G(batch_x, batch_y, seed)
        mse = F.mse_loss(x, x_hat)
        prior_loss = KLD_loss(z_mean, z_logvar)
        margin_loss = \
            kl_conditional_and_marg(z_mean, z_logvar, args.latent_dim)
        real_labels = torch.ones((batch_x.shape[0]), device=device)
        real_pred = net_D(batch_x, args.coff)
        real_loss = discriminator_loss(real_pred, real_labels)
        fake_labels = torch.zeros((x_hat.shape[0]), device=device)
        fake_pred = net_D(x_hat, args.coff)
        fake_loss = discriminator_loss(fake_pred, fake_labels)
        d_fake_loss = discriminator_loss(fake_pred, real_labels)
        # loss for Generator
        G_loss = 1 * mse + args.alpha * prior_loss + \
            args.lambda_ * (margin_loss + mse) + args.gamma * d_fake_loss
        G_loss.backward()
        optimizer_G.step()

    return models, optimizers


def train(args):
    """
    Wrapper function for training cVAE model
    """
    if args.isSaving:
        # create a folder for saving model checkpoint
        create_folder(args.checkpoint_path)
    hid_dims = []
    hid_dims.append(args.h1)
    hid_dims.append(args.h2)
    hid_dims.append(args.h3)
    hid_dims.append(args.h4)
    args.VAE_hidden_dims = hid_dims[:args.nb_layers]
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
    device = torch.device('cuda')
    # read data
    # load features we are using, here 108 ROI features
    features = txt2list(global_config.ROI_features_path)
    # read train and val data
    train_pkl = load_pkl(
        os.path.join(args.data_path, 'unmatch2match_train.pkl'))
    train_array = train_pkl['data'][features].values
    train_site_onehot = one_hot(train_pkl['data']['SITE'].values)
    train_dx_onehot = one_hot(train_pkl['data']['DX'].values, nb_classes=3)
    train_sex_onehot = one_hot(train_pkl['data']['SEX'].values - 1,
                               nb_classes=2)
    train_age = train_pkl['data'][['AGE']].values
    train_mmse = train_pkl['data'][['MMSE']].values
    train_array = torch.tensor(
        np.concatenate((train_array, train_site_onehot, train_age,
                        train_sex_onehot, train_dx_onehot, train_mmse),
                       axis=1))
    # val
    # The reason I used val_gcVAE is only using true MMSE and Diagnosis
    val_pkl = load_pkl(os.path.join(args.data_path, 'unmatch2match_val.pkl'))
    val_array = val_pkl['data'][features].values
    # !! For validation, still map with correct site
    val_site_onehot = one_hot(val_pkl['data']['SITE'].values)
    val_dx_onehot = one_hot(val_pkl['data']['DX'].values, nb_classes=3)
    val_sex_onehot = one_hot(val_pkl['data']['SEX'].values - 1, nb_classes=2)
    val_age = val_pkl['data'][['AGE']].values
    val_mmse = val_pkl['data'][['MMSE']].values
    val_array = torch.tensor(
        np.concatenate((val_array, val_site_onehot, val_age, val_sex_onehot,
                        val_dx_onehot, val_mmse),
                       axis=1)).float()
    val_array = val_array.to(device)
    # create dataloader
    dataloader = train_dataloader(train_array, args.batch_size)
    # build model, Generator(VAE) and Discriminator (D)
    net_G = VAE(in_dim=args.in_dim,
                nb_classes=args.nb_classes,
                latent_dim=args.latent_dim,
                p_dropout=args.drop_out,
                hidden_dims=args.VAE_hidden_dims)
    net_D = Discriminator(
        in_dim=args.in_dim,
        nb_classes=2,  # x or x_hat
        p_dropout=args.drop_out,
        hidden_dims=args.D_hidden_dims)

    # move to device
    net_G.to(device)
    net_D.to(device)
    models = (net_G, net_D)
    # optimizers for G and D
    # I used optimizer default settings (i.e. betas, eps) of Pytorch
    optimizer_G = torch.optim.Adam(params=net_G.parameters(),
                                   lr=args.lr,
                                   betas=(0.9, 0.999),
                                   eps=1e-7,
                                   amsgrad=False)
    optimizer_D = torch.optim.Adam(params=net_D.parameters(),
                                   lr=args.dis_lr,
                                   betas=(0.9, 0.999),
                                   eps=1e-7,
                                   amsgrad=False)
    optimizers = (optimizer_G, optimizer_D)
    scheduler_G = torch.optim.lr_scheduler.MultiStepLR(
        optimizer_G, milestones=[args.lr_step], gamma=0.1)
    scheduler_D = torch.optim.lr_scheduler.MultiStepLR(optimizer_D,
                                                       milestones=[100, 800],
                                                       gamma=0.1)
    # val mse logger
    best_G = None
    best_D = None
    best_val_mse = 1e5
    val_mse_logger = []
    # begin training
    for epoch in range(args.epochs):
        args.step = epoch
        net_G.train()
        net_D.train()
        # train 1 epoch
        models, optimizers = \
            train_1epoch(args, dataloader, models, optimizers, device)
        scheduler_G.step()
        scheduler_D.step()
        net_G.eval()
        net_D.eval()
    # since cVAE has a sampling process => which means it is random
    # performance on validation set needs to be average 100 times
    valROIs_hat_mean = torch.zeros_like(val_array[:, :args.in_dim])
    for i in range(100):
        [val_x_hat, val_x, _, _] = net_G(val_array[:, :args.in_dim],
                                         val_array[:, args.in_dim:], i)
        valROIs_hat_mean += val_x_hat / 100

    val_mse = F.mse_loss(valROIs_hat_mean, val_array[:, :args.in_dim])
    val_mse = val_mse.data.cpu().numpy()
    best_val_mse = val_mse
    val_mse_logger.append(val_mse)
    best_G = net_G
    best_D = net_D
    # save vae_mse_logger
    if args.isSaving:
        model_ckpt = args.model_name + '.pt'
        discriminator_ckpt = args.model_name + '_discrinminator.pt'
        torch.save(best_G, os.path.join(args.checkpoint_path, model_ckpt))
        torch.save(best_D,
                   os.path.join(args.checkpoint_path, discriminator_ckpt))
        list2txt(val_mse_logger,
                 os.path.join(args.checkpoint_path, 'val_mse_log.txt'))
    else:
        return best_G, best_val_mse


if __name__ == '__main__':
    train(train_covae_args_parser())
