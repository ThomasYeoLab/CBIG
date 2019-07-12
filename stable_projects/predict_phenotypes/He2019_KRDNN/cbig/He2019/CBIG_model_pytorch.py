#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch
import torch.utils.data
import torch.nn as nn
import torch.nn.init as init


class CBIG_dataset(torch.utils.data.Dataset):
    """PyTorch dataset class

    Attributes:
        x (torch.tensor): tensor for x data
        y (torch.tensor): tensor for y data
    """

    def __init__(self, x, y, for_sex=False):
        """initialization of PyTorch dataset class

        Args:
            x (ndarray): x data
            y (ndarray): y data
            for_sex (bool, optional): whether the network is used for sex
                prediction
        """
        self.x = torch.from_numpy(x).float()
        if for_sex:
            self.y = torch.squeeze(torch.from_numpy(y).long())
        else:
            self.y = torch.from_numpy(y).float().view(-1, 1)

    def __getitem__(self, idx):
        x = self.x[idx]
        y = self.y[idx]
        return x, y

    def __len__(self):
        return int(self.x.shape[0])


class fnn_3l(nn.Module):
    '''3 layer FNN model

    Attributes:
        fc1 (torch.nn.Sequential): First layer of FNN
        fc2 (torch.nn.Sequential): Second layer of FNN
        fc3 (torch.nn.Sequential): Third layer of FNN
    '''

    def __init__(self, input_size, n_l1, n_l2, dropout, for_sex=False):
        """initialization of 3 layer FNN

        Args:
            input_size (int): dimension of input data
            n_l1 (int): number of node in first layer
            n_l2 (int): number of node in second layer
            dropout (float): rate of dropout
            for_sex (bool, optional): whether the network is used for sex
                prediction
        """
        super(fnn_3l, self).__init__()
        self.fc1 = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_size, n_l1),
            nn.ReLU(),
            nn.BatchNorm1d(n_l1),
        )
        self.fc2 = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(n_l1, n_l2),
            nn.ReLU(),
            nn.BatchNorm1d(n_l2),
        )
        if for_sex:
            self.fc3 = nn.Sequential(
                nn.Dropout(p=dropout),
                nn.Linear(n_l2, 2),
            )
        else:
            self.fc3 = nn.Sequential(
                nn.Dropout(p=dropout),
                nn.Linear(n_l2, 1),
            )
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.Conv1d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                init.xavier_uniform_(m.weight)

    def forward(self, x):
        x = self.fc1(x)
        x = self.fc2(x)
        x = self.fc3(x)
        return x


class fnn_2l(nn.Module):
    '''2 layer FNN model

    Attributes:
        fc1 (torch.nn.Sequential): First layer of FNN
        fc2 (torch.nn.Sequential): Second layer of FNN
    '''

    def __init__(self, input_size, n_l1, dropout, for_sex=False):
        """initialization of 2 layer FNN

        Args:
            input_size (int): dimension of input data
            n_l1 (int): number of node in first layer
            dropout (float): rate of dropout
            for_sex (bool, optional): whether the network is used for sex
                prediction
        """
        super(fnn_2l, self).__init__()
        self.fc1 = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_size, n_l1),
            nn.ReLU(),
            nn.BatchNorm1d(n_l1),
        )
        if for_sex:
            self.fc2 = nn.Sequential(
                nn.Dropout(p=dropout),
                nn.Linear(n_l1, 2),
            )
        else:
            self.fc2 = nn.Sequential(
                nn.Dropout(p=dropout),
                nn.Linear(n_l1, 1),
            )
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.Conv1d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                init.xavier_uniform_(m.weight)

    def forward(self, x):
        x = self.fc1(x)
        x = self.fc2(x)
        return x


class Edge2Edge(nn.Module):
    """BrainNetCNN edge to edge (e2e) layer

    Attributes:
        channel (int): number of input channel
        col_conv (nn.Conv2d): column convolution
        dim (int): number of ROI for functional connectivity
        filters (int): number of output channel
        row_conv ((nn.Conv2d): row convolution
    """

    def __init__(self, channel, dim, filters):
        """initialization function of e2e layer

        Args:
            channel (int): number of input channel
            dim (int): number of ROI for functional connectivity
            filters (int): number of output channel
        """
        super(Edge2Edge, self).__init__()
        self.channel = channel
        self.dim = dim
        self.filters = filters
        self.row_conv = nn.Conv2d(channel, filters, (1, dim))
        self.col_conv = nn.Conv2d(channel, filters, (dim, 1))

    def forward(self, x):
        """e2e by two conv2d with line filter

        Args:
            x (tensor): input tensor

        Returns:
            tensor: e2e output
        """
        size = x.size()
        row = self.row_conv(x)
        col = self.col_conv(x)
        row_ex = row.expand(size[0], self.filters, self.dim, self.dim)
        col_ex = col.expand(size[0], self.filters, self.dim, self.dim)
        return row_ex + col_ex


class Edge2Node(nn.Module):
    """BrainNetCNN edge to node (e2n) layer

    Attributes:
        channel (int): number of input channel
        col_conv (nn.Conv2d): column convolution
        dim (int): number of ROI for functional connectivity
        filters (int): number of output channel
        row_conv ((nn.Conv2d): row convolution
    """

    def __init__(self, channel, dim, filters):
        """initialization function of e2n layer

        Args:
            channel (int): number of input channel
            dim (int): number of ROI for functional connectivity
            filters (int): number of output channel
        """
        super(Edge2Node, self).__init__()
        self.channel = channel
        self.dim = dim
        self.filters = filters
        self.row_conv = nn.Conv2d(channel, filters, (1, dim))
        self.col_conv = nn.Conv2d(channel, filters, (dim, 1))

    def forward(self, x):
        """e2n by add two conv2d

        Args:
            x (tensor): input tensor

        Returns:
            tensor: e2n output
        """
        row = self.row_conv(x)
        col = self.col_conv(x)
        return row + col.permute(0, 1, 3, 2)


class Node2Graph(nn.Module):
    """BrainNetCNN node to graph (n2g) layer

    Attributes:
        channel (int): number of input channel
        conv (nn.Conv2d): convolution for n2g
        dim (int): number of ROI for functional connectivity
        filters (int): number of output channel
    """

    def __init__(self, channel, dim, filters):
        """initialization function of n2g layer

        Args:
            channel (int): number of input channel
            dim (int): number of ROI for functional connectivity
            filters (int): number of output channel
        """
        super(Node2Graph, self).__init__()
        self.channel = channel
        self.dim = dim
        self.filters = filters
        self.conv = nn.Conv2d(channel, filters, (dim, 1))

    def forward(self, x):
        """n2g by convolution

        Args:
            x (tensor): input tensor

        Returns:
            tensor: n2g output
        """
        return self.conv(x)


class bcn(nn.Module):
    """BrainNetCNN (bcn) network

    Attributes:
        BatchNorm (torch.nn): batch normalization layer
        dropout (torch.nn): dropout layer
        e2e (torch.nn): e2e layer
        e2n (torch.nn): e2n layer
        fc (torch.nn): fully connected layer
        n2g (torch.nn): n2g layer
        n2g_filter (int): number of n2g filter
    """

    def __init__(self, e2e, e2n, n2g, dropout, dim, for_sex=False):
        """initialization of bcn network

        Args:
            e2e (int): number of e2e filters
            e2n (int): number of e2n filters
            n2g (int): number of n2g filters
            dropout (float): dropout rate
            dim (int): dimension of functional connectivity
            for_sex (bool, optional): whether sex is the behavioral measures
        """
        super(bcn, self).__init__()
        self.n2g_filter = n2g
        self.e2e = Edge2Edge(1, dim, e2e)
        self.e2n = Edge2Node(e2e, dim, e2n)
        self.dropout = nn.Dropout(p=dropout)
        self.n2g = Node2Graph(e2n, dim, n2g)
        if for_sex:
            self.fc = nn.Linear(n2g, 2)
        else:
            self.fc = nn.Linear(n2g, 1)
        self.BatchNorm = nn.BatchNorm1d(n2g)
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.Conv1d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                init.xavier_uniform_(m.weight)

    def forward(self, x):
        x = self.e2e(x)
        x = self.dropout(x)
        x = self.e2n(x)
        x = self.dropout(x)
        x = self.n2g(x)
        x = self.dropout(x)
        x = x.view(-1, self.n2g_filter)
        x = self.fc(self.BatchNorm(x))
        return x


class bcn_tmp(nn.Module):
    """BrainNetCNN (bcn) network

    Attributes:
        BatchNorm (torch.nn): batch normalization layer
        dropout (torch.nn): dropout layer
        e2e (torch.nn): e2e layer
        e2n (torch.nn): e2n layer
        fc (torch.nn): fully connected layer
        n2g (torch.nn): n2g layer
        n2g_filter (int): number of n2g filter
    """

    def __init__(self, e2e, e2n, n2g, dropout, dim, alpha, for_sex=False):
        """initialization of bcn network

        Args:
            e2e (int): number of e2e filters
            e2n (int): number of e2n filters
            n2g (int): number of n2g filters
            dropout (float): dropout rate
            dim (int): dimension of functional connectivity
            for_sex (bool, optional): whether sex is the behavioral measures
        """
        super(bcn_tmp, self).__init__()
        self.n2g_filter = n2g
        self.e2e = Edge2Edge(1, dim, e2e)
        self.e2n = Edge2Node(e2e, dim, e2n)
        self.dropout = nn.Dropout(p=dropout)
        self.LeakyReLU = nn.LeakyReLU(alpha)
        self.n2g = Node2Graph(e2n, dim, n2g)
        if for_sex:
            self.fc = nn.Linear(n2g, 2)
        else:
            self.fc = nn.Linear(n2g, 1)
        self.BatchNorm = nn.BatchNorm1d(n2g)
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.Conv1d):
                init.xavier_uniform_(m.weight)
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                init.xavier_uniform_(m.weight)

    def forward(self, x):
        x = self.e2e(x)
        x = self.e2n(x)
        x = self.dropout(x)
        x = self.LeakyReLU(x)
        x = self.n2g(x)
        x = self.LeakyReLU(x)
        x = x.view(-1, self.n2g_filter)
        x = self.fc(self.BatchNorm(x))
        return x


def main():
    pass


if __name__ == '__main__':
    main()
