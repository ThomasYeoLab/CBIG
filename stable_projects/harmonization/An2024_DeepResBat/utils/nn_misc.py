#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import torch
from torch.utils.data import DataLoader, WeightedRandomSampler
import torch.nn as nn
from torch.nn import init


def train_dataloader(data, batch_size):
    """
    Training data loader
    Args:
        data (array): tensor array
        batch_size (int): Batch size
    """

    return DataLoader(data,
                      batch_size=batch_size,
                      shuffle=True,
                      drop_last=False,
                      num_workers=1)


def train_dataloader_balanced(data, batch_size, class_digit=108):
    """
    A balanced dataloader for two classes data
    The <class_digit> == 1 means ADNI data
    The <class_digit> == 0 means nonADNI data
    Give balanced ratio batch by oversampling

    Args:
        data (_type_): _description_
        batch_size (_type_): _description_
        class_digit (int, optional): _description_. Defaults to 108.
    """
    class_vec = data[:, class_digit].int()

    # define weights for balanced sampling

    def make_weights_for_balanced_classes(data, class_digit):
        nb_classes = len(torch.unique(data[:, class_digit]))
        count = [0] * nb_classes
        for i in range(data.shape[0]):
            count[class_vec[i]] += 1
        weight_per_class = [0.] * nb_classes
        N = float(sum(count))
        for i in range(nb_classes):
            weight_per_class[i] = N / float(count[i])
        sum_weight = sum(weight_per_class)
        for i in range(nb_classes):
            weight_per_class[i] /= sum_weight
        weight = [0] * data.shape[0]
        for i in range(data.shape[0]):
            weight[i] = weight_per_class[class_vec[i]]
        return weight

    weights = make_weights_for_balanced_classes(data, class_digit)
    weights = torch.Tensor(weights)
    sampler = WeightedRandomSampler(weights, len(weights), replacement=True)

    return DataLoader(data,
                      batch_size=batch_size,
                      shuffle=False,
                      drop_last=False,
                      sampler=sampler,
                      num_workers=1)


def xavier_uniform_init(m, init_gain=1.0):
    """
    xavier uniform initialization
    Args:
        m (class Layer): Layer
        init_gain (float, optional): Noise level
    """
    if type(m) == nn.Linear:
        init.xavier_uniform_(m.weight, gain=init_gain)
        m.bias.data.fill_(0)


def KLD_loss(mu, logvar):
    """
    KL divergence loss
    Args:
        mu (array): mean
        logvar (array): log of variance
    """

    kld_loss = -0.5 * torch.sum(1 + logvar - mu**2 - logvar.exp(), dim=1)
    kld_loss = torch.mean(kld_loss)

    return kld_loss


def all_pairs_gaussian_kl(mu, sigma, add_third_term=False):
    """
    Pair-wise Gaussian KL divergence loss
    Args:
        mu (array): mean
        sigma (array): sigma
        add_third_term (bool, optional): Whether to add third term in loss.
    """
    sigma_sq = sigma**2 + 1e-8
    sigma_sq_inv = torch.reciprocal(sigma_sq)

    # dot product of all sigma_inv vectors with sigma is the same as a matrix
    # mult of diag
    first_term = torch.matmul(sigma_sq, torch.transpose(sigma_sq_inv, 0, 1))
    r = torch.matmul(mu * mu, torch.transpose(sigma_sq_inv, 0, 1))
    r2 = mu * mu * sigma_sq_inv
    r2 = torch.sum(r2, dim=1)
    r2 = torch.reshape(r2, [-1, 1])

    # square distance
    second_term = \
        2 * torch.matmul(mu, torch.transpose(mu * sigma_sq_inv, 0, 1))
    second_term = r - second_term + torch.transpose(r2, 0, 1)

    if add_third_term:
        r = torch.sum(torch.log(sigma_sq), dim=1)
        r = torch.reshape(r, [-1, 1])
        third_term = r - torch.transpose(r, 0, 1)
    else:
        third_term = 0
    return 0.5 * (first_term + second_term + third_term)


def kl_conditional_and_marg(z_mean, z_log_sigma_sq, dim_z):
    """
    Conditional and marginal KL divergence loss
    Args:
        z_mean (array): mean of z
        z_log_sigma_sq (array): Log of squaared sigma of z
        dim_z (int): Dimension of z
    """
    z_sigma = torch.exp(0.5 * z_log_sigma_sq)
    all_pairs_GKL = all_pairs_gaussian_kl(z_mean, z_sigma, False) - 0.5 * dim_z
    return torch.mean(all_pairs_GKL)


def discriminator_loss(pred, labels):
    """
    Loss for Discriminator
        1. pred is a probability matrix after LogSigmoid() function
        Shape of pred is[batch_size, 2], 2 is for real or fake
        2. labels is a one-hot vector, length is batch_size
    Args:
        pred (array): Prediction
        labels (array): Ground truth
    """
    dis_loss = torch.nn.NLLLoss(reduction='mean')
    labels = labels.long()

    return dis_loss(pred, labels)


def vae_harm_predict(model, x, seed, site=[0, 0], recon=False):
    """
    Mapping for VAE models

    Args:
        model (_type_): _description_
        data (_type_): _description_
        seed (_type_): _description_
        site (list, optional): _description_. Defaults to [0, 0].
    """
    sites = x[:, 108:110]
    rois = x[:, :108]
    if recon:
        pass
    else:
        sites[:, 0] = site[0]
        sites[:, 1] = site[1]
    covs = x[:, 110:]
    all_covs = torch.cat((sites, covs), dim=1)
    pred = (model(rois, all_covs, seed)[0]).detach().numpy()

    return pred
