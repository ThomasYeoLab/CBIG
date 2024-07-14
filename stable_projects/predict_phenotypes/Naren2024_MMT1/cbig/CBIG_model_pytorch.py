#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch
import torch.nn as nn
import numpy as np
import nibabel as nib
from cbig.config import config

if torch.cuda.is_available():
    data_type = torch.cuda.FloatTensor
else:
    data_type = torch.FloatTensor


def torch_nanmean(x, mask):
    '''Calculate mean and omit NAN

    Args:
        x (torch.tensor): input data
        mask (torch.tensor, optional): mask indicated NAN

    Returns:
        torch.Tensor: mean value (omit NAN)
    '''

    num = torch.where(mask, torch.full_like(x, 0), torch.full_like(x, 1)).sum()
    value = torch.where(mask, torch.full_like(x, 0), x).sum()
    return value / num


def msenanloss(input, target, mask=None):
    '''Calculate MSE (mean absolute error) and omit NAN

    Args:
        input (torch.tensor): predicted value
        target (torch.tensor): original value
        mask (torch.tensor, optional): mask indicated NAN

    Returns:
        torch.Tensor: MSE loss (omit NAN)
    '''

    ret = (input - target)**2
    if mask is None:
        mask = torch.isnan(ret)
    return torch_nanmean(ret, mask)


def crop_center(data, out_sp):
    '''Crop 3D volumetric input size

    Args:
        data (ndarray): input data (182, 218, 182)
        out_sp (tuple): output size (160, 192, 160)

    Returns:
        data_crop (ndarray): cropped data
    '''

    in_sp = data.shape
    nd = np.ndim(data)
    x_crop = int((in_sp[-1] - out_sp[-1]) / 2)
    y_crop = int((in_sp[-2] - out_sp[-2]) / 2)
    z_crop = int((in_sp[-3] - out_sp[-3]) / 2)
    if nd == 3:
        data_crop = data[x_crop:-x_crop, y_crop:-y_crop, z_crop:-z_crop]
    elif nd == 4:
        data_crop = data[:, x_crop:-x_crop, y_crop:-y_crop, z_crop:-z_crop]
    else:
        raise ('Wrong dimension! dim=%d.' % nd)
    return data_crop


class vol_dataset(torch.utils.data.Dataset):
    """PyTorch dataset class for volumetric data

    Attributes:
        x (ndarray): volumetric data
        y (ndarray): phenotype data
        icv (ndarray): intracerebroventricular (ICV)
    """

    def __init__(self, x, y, icv=None):
        """initialization of PyTorch dataset class
        """
        self.sublist = x
        self.y = torch.from_numpy(y).float()
        self.cutoff_size = config.CUTOFF_SIZE
        self.icv = icv

        if icv is not None:
            self.icv = torch.from_numpy(icv).float()[:, None, None, None]

    def __getitem__(self, idx):
        sublist = self.sublist[idx]
        nifti_img = nib.load(sublist)
        nii_data = nifti_img.get_fdata()
        # Preprocessing
        nii_data = crop_center(nii_data, self.cutoff_size)
        nii_data[nii_data < 0] = 0
        x = torch.unsqueeze(
            torch.from_numpy(nii_data / nii_data.mean()).float(), 0)

        y = self.y[idx]

        if self.icv is not None:
            icv = self.icv[idx]
            return x, y, icv
        else:
            return x, y

    def __len__(self):
        return int(self.sublist.shape[0])


class SFCN(nn.Module):
    ''' Simple Fully Convolutional Network

    Attributes:
        channel_number (list): channel number of each convolution layer
        output_dim (int): output dimensionality of SFCN
        dropout (float): dropout rate
        feature_extractor (torch.nn.Sequential): feature extractior of SFCN
        classifier (torch.nn.Sequential): classifier of SFCN
    '''

    def __init__(self,
                 channel_number=config.CHANNEL_NUMBER,
                 output_dim=config.OUTPUT_DIM,
                 dropout=None):
        super(SFCN, self).__init__()

        n_layer = len(channel_number)
        self.feature_extractor = nn.Sequential()
        for i in range(n_layer):
            if i == 0:
                in_channel = 1
            else:
                in_channel = channel_number[i - 1]
            out_channel = channel_number[i]
            if i < n_layer - 1:
                self.feature_extractor.add_module(
                    'conv_%d' % i,
                    self.conv_layer(in_channel,
                                    out_channel,
                                    maxpool=True,
                                    kernel_size=3,
                                    padding=1))
            else:
                self.feature_extractor.add_module(
                    'conv_%d' % i,
                    self.conv_layer(in_channel,
                                    out_channel,
                                    maxpool=False,
                                    kernel_size=1,
                                    padding=0))
        self.classifier = nn.Sequential()
        avg_shape = [5, 6, 5]
        self.classifier.add_module('average_pool', nn.AvgPool3d(avg_shape))
        if dropout:
            self.classifier.add_module('dropout', nn.Dropout(dropout))
        i = n_layer
        in_channel = channel_number[-1] + 1
        out_channel = output_dim
        self.classifier.add_module(
            'conv_%d' % i,
            nn.Conv3d(in_channel, out_channel, padding=0, kernel_size=1))

    @staticmethod
    def conv_layer(in_channel,
                   out_channel,
                   maxpool=True,
                   kernel_size=3,
                   padding=0,
                   maxpool_stride=2):
        if maxpool is True:
            layer = nn.Sequential(
                nn.Conv3d(in_channel,
                          out_channel,
                          padding=padding,
                          kernel_size=kernel_size),
                nn.BatchNorm3d(out_channel),
                nn.MaxPool3d(2, stride=maxpool_stride),
                nn.ReLU(),
            )
        else:
            layer = nn.Sequential(
                nn.Conv3d(in_channel,
                          out_channel,
                          padding=padding,
                          kernel_size=kernel_size),
                nn.BatchNorm3d(out_channel), nn.ReLU())
        return layer

    def forward(self, x, icv):

        x = self.feature_extractor(x)
        x = self.classifier.average_pool(x)
        x = torch.cat((x, icv), 1)
        x = self.classifier.dropout(x)
        x = self.classifier.conv_6(x)

        return x.reshape(x.shape[0], -1)


class transfer_model(nn.Module):
    ''' Model for finetuning classical transfer learning or
        meta-matching finetune (for speedup training)

    Attributes:
        pretrained_net (torch.nn.Sequential): model of pre-trained SFCN
        output_dim (int): output dimensionality of the transfer model
        dropout (float): dropout rate
        initial_last_layer (bool): whether initialize last layer
        with parameters from pre-trained model
        param (float): selected parameters from pre-trained model
        in_c (int): input channel number
        out_c (int): output channel number
    '''

    def __init__(self,
                 pretrained_net,
                 output_dim=config.OUTPUT_DIM_TRANS,
                 dropout=config.DROPOUT,
                 initial_last_layer=False,
                 param=None,
                 in_c=config.IN_C_TRANS,
                 out_c=config.OUT_C_TRANS):
        super(transfer_model, self).__init__()

        in_channel = in_c
        out_channel = out_c
        layer_idx = config.LAYER_IDX_TRANS
        self.feature_extractor = nn.Sequential()
        self.feature_extractor.add_module(
            'conv_%d' % layer_idx,
            nn.Sequential(
                nn.Conv3d(in_channel, out_channel, padding=0, kernel_size=1),
                nn.BatchNorm3d(out_channel), nn.ReLU()))

        self.classifier = nn.Sequential()
        avg_shape = [5, 6, 5]
        self.classifier.add_module('average_pool', nn.AvgPool3d(avg_shape))
        if dropout:
            self.classifier.add_module('dropout', nn.Dropout(dropout))
        i = config.FLAYER_IDX_TRANS
        in_channel = int(out_c + 1)
        out_channel = output_dim
        self.classifier.add_module(
            'conv_%d' % i,
            nn.Conv3d(in_channel, out_channel, padding=0, kernel_size=1))

        self.initialize(pretrained_net, initial_last_layer, param)

    def initialize(self, pretrained_net, initial_last_layer, param):

        penultimate_layer = pretrained_net.feature_extractor.conv_5
        tmp_weight = penultimate_layer[0].weight
        tmp_bias = penultimate_layer[0].bias
        self.feature_extractor.conv_5[0].weight = nn.Parameter(tmp_weight)
        self.feature_extractor.conv_5[0].bias = nn.Parameter(tmp_bias)

        if initial_last_layer:
            fc_last = pretrained_net.classifier.conv_6
            tmp_weight = fc_last.weight[param, :].unsqueeze(0)
            tmp_bias = fc_last.bias[param].unsqueeze(0)
            self.classifier.conv_6.weight = nn.Parameter(tmp_weight)
            self.classifier.conv_6.bias = nn.Parameter(tmp_bias)

    def forward(self, x):
        x, icv = x[:, :-1, :, :, :].type(data_type), torch.reshape(
            x[:, -1, :, :, :],
            (x.shape[0], int(
                x.shape[2] * x.shape[3] * x.shape[4])))[:, 0].type(data_type)

        x = self.feature_extractor(x)
        icv = icv[:, None, None, None, None]
        x = self.classifier.average_pool(x)
        x = torch.cat((x, icv), 1)
        x = self.classifier.dropout(x)
        x = self.classifier.conv_6(x)

        return x.reshape(x.shape[0], -1)
