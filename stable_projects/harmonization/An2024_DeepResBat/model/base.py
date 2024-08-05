#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
from torch import nn
from abc import abstractmethod
from typing import TypeVar


class BaseFCN(nn.Module):
    """
    Model base for FCN (fully connected neural network)
    """

    def __init__(self):
        super(BaseFCN, self).__init__()

    @abstractmethod
    def forward(self, input):
        pass

    @abstractmethod
    def init_parameters(self):
        pass


Tensor = TypeVar('torch.tensor')


class BaseXGB:

    def __init__(self) -> None:
        pass
