#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
from torch import nn
from model.base import BaseFCN
from utils.nn_misc import xavier_uniform_init


class goalDNN(BaseFCN):
    """
    goalDNN model
    """

    def __init__(self, in_dim, nb_category, nb_measures, p_dropout,
                 hidden_dims):
        super(goalDNN, self).__init__()
        modules = []
        if hidden_dims is None:
            hidden_dims = [512, 256]
        # build model
        for h_dim in hidden_dims:
            modules.append(
                nn.Sequential(
                    nn.Dropout(p=p_dropout), nn.Linear(in_dim, h_dim),
                    nn.ReLU()))
            in_dim = h_dim
        self.model = nn.Sequential(*modules)
        # output layers
        self.hid2cat = nn.Linear(hidden_dims[-1], nb_category)
        self.hid2measure = nn.Linear(hidden_dims[-1], nb_measures)

        # initialization
        self.init_parameters()

    def forward(self, input):
        """
        Forward pass of TaskDNN
        """
        hid = self.model(input)
        cat = self.hid2cat(hid)
        measure = self.hid2measure(hid)

        return [measure, cat]

    def init_parameters(self):
        """
        weight and bias initialization
        """
        self.model.apply(xavier_uniform_init)
        self.hid2measure.apply(xavier_uniform_init)
        self.hid2cat.apply(xavier_uniform_init)
