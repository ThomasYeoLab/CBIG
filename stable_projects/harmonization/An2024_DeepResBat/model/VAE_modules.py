#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import torch
from torch import nn
from model.base import BaseFCN, Tensor
from typing import List
from utils.nn_misc import xavier_uniform_init


class Encoder(BaseFCN):

    def __init__(self,
                 in_dim,
                 latent_dim,
                 p_dropout,
                 hidden_dims: List = None) -> None:
        super(Encoder, self).__init__()

        modules = []
        if hidden_dims is None:
            hidden_dims = [256, 128]

        # Build Encoder
        for h_dim in hidden_dims:
            modules.append(
                nn.Sequential(nn.Dropout(p=p_dropout),
                              nn.Linear(in_dim, h_dim), nn.Tanh()))
            in_dim = h_dim

        self.encoder = nn.Sequential(*modules)
        self.fc_mu = nn.Linear(hidden_dims[-1], latent_dim)
        self.fc_logvar = nn.Linear(hidden_dims[-1], latent_dim)
        # initialization
        self.init_parameters()

    def forward(self, input: Tensor) -> List[Tensor]:
        """
        Encoder
        """
        result = self.encoder(input)
        mu = self.fc_mu(result)
        mu = torch.tanh(mu)
        logvar = self.fc_logvar(result)

        return [mu, logvar]

    def init_parameters(self):
        """
        weight and bias initialization
        """
        self.encoder.apply(xavier_uniform_init)
        self.fc_mu.apply(xavier_uniform_init)
        self.fc_logvar.apply(xavier_uniform_init)


class Decoder(BaseFCN):

    def __init__(self,
                 in_dim,
                 nb_clases,
                 latent_dim,
                 p_dropout,
                 hidden_dims: List = None) -> None:
        super(Decoder, self).__init__()

        modules = []
        if hidden_dims is None:
            hidden_dims = [256, 128]

        self.decoder_input = nn.Linear(latent_dim + nb_clases, hidden_dims[-1])
        hidden_dims.reverse()

        # Build Decoder
        for i in range(len(hidden_dims) - 1):
            modules.append(
                nn.Sequential(nn.Dropout(p=p_dropout),
                              nn.Linear(hidden_dims[i], hidden_dims[i + 1]),
                              nn.Tanh()))

        self.decoder = nn.Sequential(*modules)
        # final layer
        self.final_layer = nn.Linear(hidden_dims[-1], in_dim)
        # initialization
        self.init_parameters()

    def forward(self, z: Tensor) -> List[Tensor]:
        """
        Decoder
        """
        result = self.decoder_input(z)
        result = torch.tanh(result)
        result = self.decoder(result)
        result = self.final_layer(result)
        return result

    def init_parameters(self):
        """
        weight and bias initialization
        """
        self.decoder_input.apply(xavier_uniform_init)
        self.decoder.apply(xavier_uniform_init)
        self.final_layer.apply(xavier_uniform_init)


class DeepComBatDecoder(BaseFCN):

    def __init__(self,
                 in_dim,
                 covs_dim,
                 nb_clases,
                 latent_dim,
                 p_dropout,
                 hidden_dims: List = None) -> None:
        super(DeepComBatDecoder, self).__init__()

        modules = []
        if hidden_dims is None:
            hidden_dims = [256, 128]

        self.decoder_input = nn.Linear(latent_dim + covs_dim + nb_clases,
                                       hidden_dims[-1])
        hidden_dims.reverse()

        # Build Decoder
        for i in range(len(hidden_dims) - 1):
            modules.append(
                nn.Sequential(nn.Dropout(p=p_dropout),
                              nn.Linear(hidden_dims[i], hidden_dims[i + 1]),
                              nn.Tanh()))

        self.decoder = nn.Sequential(*modules)
        # final layer
        self.final_layer = nn.Linear(hidden_dims[-1], in_dim)
        # initialization
        self.init_parameters()

    def forward(self, z: Tensor) -> List[Tensor]:
        """
        Decoder
        """
        result = self.decoder_input(z)
        result = torch.tanh(result)
        result = self.decoder(result)
        result = self.final_layer(result)
        return result

    def init_parameters(self):
        """
        weight and bias initialization
        """
        self.decoder_input.apply(xavier_uniform_init)
        self.decoder.apply(xavier_uniform_init)
        self.final_layer.apply(xavier_uniform_init)


class Discriminator(BaseFCN):

    def __init__(self,
                 in_dim,
                 nb_classes,
                 p_dropout,
                 hidden_dims: List = None) -> None:
        super(Discriminator, self).__init__()

        modules = []
        if hidden_dims is None:
            hidden_dims = [32, 32]

        # Build Discriminator
        for h_dim in hidden_dims:
            modules.append(
                nn.Sequential(nn.Dropout(p=p_dropout),
                              nn.Linear(in_dim, h_dim), nn.Tanh()))
            in_dim = h_dim
        self.discriminator = nn.Sequential(*modules)
        # final layer
        self.final_layer = nn.Linear(hidden_dims[-1], nb_classes)
        # initialization
        self.init_parameters()

    def forward(self, x: Tensor, coff) -> Tensor:
        """
        Forward pass
        """
        # Get the prediction of input, real or fake
        x = x + coff * torch.randn_like(x)
        pred = self.discriminator(x)
        pred = self.final_layer(pred)
        # binary classification, distinguish input is x or x_hat
        pred = torch.log_softmax(pred, dim=1)

        return pred

    def init_parameters(self):
        """
        weight and bias initialization
        """
        self.discriminator.apply(xavier_uniform_init)
        self.final_layer.apply(xavier_uniform_init)
