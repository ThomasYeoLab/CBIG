#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

from keras.layers import Input, Dropout, Dense, Activation, LeakyReLU, Conv2D
from keras.layers import Lambda, Flatten
from keras.models import Model
from keras.layers.merge import add
from keras.layers.core import Permute
from keras import backend as K
from keras.regularizers import l2
from keras.layers.normalization import BatchNormalization
from graph import GraphConvolution


def gcnn(input_shape,
         dropout,
         n_l1,
         graph_matrix,
         support,
         l2_reg,
         for_sex=False,
         n_measure=None):
    """GCNN model

    Args:
        input_shape (int): dimension of input x
        dropout (float): dropout rate
        n_l1 (int): number of node for GCNN layer 1
        graph_matrix (ndarray): graph matrix for GCNN calculation
        support (int): support for filter
        l2_reg (float): l2 regularizer rate
        for_sex (bool, optional): whether the network is used for sex
            prediction
        n_measure (int, optional): number of measures (HCP dataset only)

    Returns:
        keras.models.Model: GCNN model
    """
    model_in = Input(shape=(input_shape, ))
    H = Dropout(dropout)(model_in)
    H = GraphConvolution(
        n_l1,
        graph_matrix=graph_matrix,
        support=support,
        activation='relu',
        W_regularizer=l2(l2_reg))(H)
    H = BatchNormalization()(H)
    H = Dropout(dropout)(H)
    if n_measure is not None:
        model_out = GraphConvolution(
            n_measure, graph_matrix=graph_matrix, support=support)(H)
    elif for_sex:
        model_out = GraphConvolution(
            2,
            graph_matrix=graph_matrix,
            support=support,
            activation='softmax')(H)
    else:
        model_out = GraphConvolution(
            1, graph_matrix=graph_matrix, support=support)(H)

    return Model(model_in, model_out)


def hcp_fnn(input_shape, n_measure, n_l1, n_l2, n_l3, dropout, l2_reg):
    """FNN model for HCP dataset

    Args:
        input_shape (int): dimension of input x
        n_measure (int): number of behavioral measures
        n_l1 (int): number of node for GCNN layer 1
        n_l2 (int): number of node for GCNN layer 2
        n_l3 (int): number of node for GCNN layer 3
        dropout (float): dropout rate
        l2_reg (float): l2 regularizer rate

    Returns:
        keras.models.Model: FNN model
    """
    init_method = 'glorot_uniform'
    model_in = Input(shape=(input_shape, ))
    H = Dropout(dropout)(model_in)
    H = Dense(n_l1, kernel_initializer=init_method)(H)
    H = Activation('linear')(H)
    H = BatchNormalization()(H)
    H = Dropout(dropout)(H)
    H = Dense(
        n_l2, kernel_initializer=init_method, kernel_regularizer=l2(l2_reg))(H)
    H = Activation('linear')(H)
    H = BatchNormalization()(H)
    H = Dropout(dropout)(H)
    H = Dense(n_l3, kernel_initializer=init_method)(H)
    H = Activation('linear')(H)
    H = BatchNormalization()(H)
    H = Dropout(dropout)(H)
    model_out = Dense(n_measure, kernel_initializer=init_method)(H)
    return Model(model_in, model_out)


def hcp_Edge2Edge(inputs, dim, filters, activation):
    """BrainNetCNN edge to edge (e2e) layer for HCP dataset

    Args:
        inputs (keras.layer): input keras layer
        dim (int): number of ROI for functional connectivity
        filters (int): number of output channel
        activation (keras.layer.Activation): activation function

    Returns:
        keras.layer: e2e keras output
    """
    row = Conv2D(filters, (1, dim), padding='valid', activation=None)(inputs)
    col = Conv2D(filters, (dim, 1), padding='valid', activation=None)(inputs)
    tiled_row = Lambda(lambda x: K.repeat_elements(x, rep=dim, axis=2),
                       (dim, dim, filters))(row)
    tiled_col = Lambda(lambda x: K.repeat_elements(x, rep=dim, axis=1),
                       (dim, dim, filters))(col)

    Sum = add([tiled_row, tiled_col])
    return activation(Sum)


def hcp_Edge2Node(inputs, dim, filters, activation):
    """BrainNetCNN edge to node (e2n) layer for HCP dataset

    Args:
        inputs (keras.layer): input keras layer
        dim (int): number of ROI for functional connectivity
        filters (int): number of output channel
        activation (keras.layer.Activation): activation function

    Returns:
        keras.layer: e2n keras output
    """
    row = Conv2D(filters, (1, dim), padding='valid', activation=None)(inputs)
    col = Conv2D(filters, (dim, 1), padding='valid', activation=None)(inputs)

    Sum = add([row, Permute(dims=(2, 1, 3))(col)])
    return activation(Sum)


def hcp_Node2Graph(inputs, dim, filters, activation):
    """BrainNetCNN node to graph (n2g) layer for HCP dataset

    Args:
        inputs (keras.layer): input keras layer
        dim (int): number of ROI for functional connectivity
        filters (int): number of output channel
        activation (keras.layer.Activation): activation function

    Returns:
        keras.layer: n2g keras output
    """
    nodes = Conv2D(filters, (dim, 1), padding='valid', activation=None)(inputs)

    return activation(nodes)


def hcp_brainnetcnn(dim,
                    n_measure,
                    e2e,
                    e2n,
                    n2g,
                    dropout,
                    leaky_alpha,
                    nb_features=1):
    """BrainNetCNN network for HCP datasets

    Args:
        dim (int): number of ROI for functional connectivity
        n_measure (int): number of behavioral measures
        e2e (int): number of e2e filters
        e2n (int): number of e2n filters
        n2g (int): number of n2g filters
        dropout (float): dropout rate
        leaky_alpha (float): leaky ReLU alpha rate
        nb_features (int, optional): number of features for input

    Returns:
        keras.models.Model: BrainNetCNN model
    """
    activation = LeakyReLU(alpha=leaky_alpha)
    In = Input(shape=(dim, dim, nb_features))
    layer_beta = hcp_Edge2Edge(In, dim, e2e, activation=activation)
    layer0 = hcp_Edge2Node(layer_beta, dim, e2n, activation=activation)
    layer1 = Dropout(dropout)(layer0)
    layer2 = hcp_Node2Graph(layer1, dim, n2g, activation=activation)
    layer3 = Flatten()(layer2)
    layer4 = Dense(n_measure, activation='linear')(layer3)
    return Model(inputs=In, outputs=layer4)


def main():
    pass


if __name__ == '__main__':
    main()
