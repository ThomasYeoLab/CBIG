#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import torch
from torch import nn
from typing import List
from model.base import Tensor
from model.VAE_modules import Encoder, Decoder, DeepComBatDecoder


class VAE(nn.Module):

    def __init__(self,
                 in_dim,
                 nb_classes,
                 latent_dim,
                 p_dropout,
                 hidden_dims: List = None) -> None:
        super(VAE, self).__init__()
        self.in_dim = in_dim
        self.nb_classes = nb_classes
        self.latent_dim = latent_dim
        self.p_dropout = p_dropout
        self.hidden_dims = hidden_dims
        self.encoder = Encoder(self.in_dim, self.latent_dim, self.p_dropout,
                               self.hidden_dims)
        self.decoder = Decoder(self.in_dim, self.nb_classes, self.latent_dim,
                               self.p_dropout, self.hidden_dims)

    def reparameterize(self, mu: Tensor, logvar: Tensor, seed: int) -> Tensor:
        """
        Reparameterize trick
        """
        std = torch.exp(0.5 * logvar)
        torch.manual_seed(seed)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x: Tensor, label: Tensor, seed: int,
                **kwargs) -> List[Tensor]:
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar, seed)
        catted_z = torch.cat([z, label], dim=1)

        recon = self.decoder(catted_z)

        return [recon, x, mu, logvar]


class DeepComBat(nn.Module):

    def __init__(self,
                 in_dim,
                 covs_dim,
                 nb_classes,
                 latent_dim,
                 p_dropout,
                 hidden_dims: List = None) -> None:
        super(DeepComBat, self).__init__()
        self.in_dim = in_dim
        self.covs_dim = covs_dim
        self.nb_classes = nb_classes
        self.latent_dim = latent_dim
        self.p_dropout = p_dropout
        self.hidden_dims = hidden_dims
        self.encoder = Encoder(self.in_dim + self.covs_dim + self.nb_classes,
                               self.latent_dim, self.p_dropout,
                               self.hidden_dims)
        self.decoder = DeepComBatDecoder(self.in_dim, self.covs_dim,
                                         self.nb_classes, self.latent_dim,
                                         self.p_dropout, self.hidden_dims)

    def reparameterize(self, mu: Tensor, logvar: Tensor, seed: int) -> Tensor:
        """
        Reparameterize trick
        """
        std = torch.exp(0.5 * logvar)
        torch.manual_seed(seed)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x: Tensor, covs: Tensor, label: Tensor, seed: int,
                **kwargs) -> List[Tensor]:
        mu, logvar = self.encoder(torch.cat([x, covs, label], dim=1))
        z = self.reparameterize(mu, logvar, seed)
        catted_z = torch.cat([z, covs, label], dim=1)

        recon = self.decoder(catted_z)

        return [recon, x, mu, logvar]
