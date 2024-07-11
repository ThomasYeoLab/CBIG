#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Pansheng Chen and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch
import torch.nn as nn
import torch.utils.data
import torch.nn.init as init


def torch_nanmean(x, mask):
    '''Calculate mean and omit NAN

    Args:
        x (torch.tensor): input data
        mask (torch.tensor): mask indicated NAN

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


class multi_task_dataset(torch.utils.data.Dataset):
    """PyTorch dataset class

    Attributes:
        x (torch.tensor): tensor for x data
        y (torch.tensor): tensor for y data
    """

    def __init__(self, x, y, for_finetune=False):
        """initialization of PyTorch dataset class

        Args:
            x (ndarray): x data
            y (ndarray): y data
            for_finetune (bool, optional): whether the network is used for
                finetune
        """
        self.x = torch.from_numpy(x).float()
        if for_finetune:
            self.y = torch.from_numpy(y).float().view(-1, 1)
        else:
            self.y = torch.from_numpy(y).float()

    def __getitem__(self, idx):
        x = self.x[idx]
        y = self.y[idx]
        return x, y

    def __len__(self):
        return int(self.x.shape[0])


class dnn(nn.Module):
    '''DNN model
    '''

    def __init__(self,
                 input_size,
                 n_hidden_layer,
                 n_l1,
                 n_l2,
                 n_l3,
                 n_l4,
                 dropout,
                 output_size=1):
        """initialization of 3 layer DNN

        Args:
            input_size (int): dimension of input data
            n_hidden_layer (int): number of hidden layers
            n_l1 (int): number of node in first layer
            n_l2 (int): number of node in second layer
            n_l3 (int): number of node in third layer
            n_l4 (int): number of node in fourth layer
            dropout (float): rate of dropout
            output_size (int, optional): dimension of output data
        """
        super(dnn, self).__init__()

        self.layers = nn.ModuleList()
        n_nodes = [n_l1, n_l2, n_l3, n_l4]
        self.layers.append(
            nn.Sequential(
                nn.Dropout(p=dropout),
                nn.Linear(input_size, n_l1),
                nn.ReLU(),
                nn.BatchNorm1d(n_l1),
            ))

        for i in range(n_hidden_layer):
            self.layers.append(
                nn.Sequential(
                    nn.Dropout(p=dropout),
                    nn.Linear(n_nodes[i], n_nodes[i + 1]),
                    nn.ReLU(),
                    nn.BatchNorm1d(n_nodes[i + 1]),
                ))
        self.layers.append(
            nn.Sequential(
                nn.Dropout(p=dropout),
                nn.Linear(n_nodes[n_hidden_layer], output_size),
            ))

        for m in self.modules():
            if isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                init.xavier_uniform_(m.weight)

    def forward(self, x):
        for layer in self.layers[:]:
            x = layer(x)

        return x
