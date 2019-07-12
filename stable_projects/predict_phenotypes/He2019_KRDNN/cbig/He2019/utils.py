#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified version of code from https://github.com/tkipf/keras-gcn
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

from __future__ import print_function

import scipy.sparse as sp
import numpy as np
from pathlib import PurePath
from scipy.sparse.linalg.eigen.arpack import eigsh, ArpackNoConvergence


def load_graph(dimension,
               path="graph/",
               graph="956Subject_corr_option_3_param_5"):
    """Load citation graph"""
    print('Loading {} graph...'.format(graph))

    # build graph
    path = PurePath(path, "{}.cites".format(graph)).as_posix()
    edges = np.genfromtxt(path, dtype=np.float32)
    index = edges[:, :2].astype(int)
    adj = sp.coo_matrix((edges[:, 2], (index[:, 0] - 1, index[:, 1] - 1)),
                        shape=(dimension, dimension),
                        dtype=np.float32)

    # build symmetric adjacency matrix
    adj = adj + adj.T.multiply(adj.T > adj) - adj.multiply(adj.T > adj)

    print('Graph has {} nodes, {} edges.'.format(adj.shape[0], edges.shape[0]))

    return adj


def normalize_adj(adj, symmetric=True):
    if symmetric:
        d = sp.diags(np.power(np.array(adj.sum(1)), -0.5).flatten(), 0)
        a_norm = adj.dot(d).transpose().dot(d).tocsr()
    else:
        d = sp.diags(np.power(np.array(adj.sum(1)), -1).flatten(), 0)
        a_norm = d.dot(adj).tocsr()
    return a_norm


def preprocess_adj(adj, symmetric=True):
    adj = adj + sp.eye(adj.shape[0])
    adj = normalize_adj(adj, symmetric)
    adj = adj.todense()
    return adj


def normalized_laplacian(adj, symmetric=True):
    adj_normalized = normalize_adj(adj, symmetric)
    laplacian = sp.eye(adj.shape[0]) - adj_normalized
    return laplacian


def rescale_laplacian(laplacian):
    try:
        print(
            'Calculating largest eigenvalue of normalized graph Laplacian...')
        largest_eigval = eigsh(
            laplacian, 1, which='LM', return_eigenvectors=False)[0]
    except ArpackNoConvergence:
        print('Eigenvalue calculation did not converge! ',
              'Using largest_eigval=2 instead.')
        largest_eigval = 2

    scaled_laplacian = (2. / largest_eigval) * laplacian - sp.eye(
        laplacian.shape[0])
    return scaled_laplacian


def chebyshev_polynomial(X, k):
    """Calculate Chebyshev polynomials up to order k.
    Return a list of sparse matrices."""
    print("Calculating Chebyshev polynomials up to order {}...".format(k))

    T_k = list()
    T_k.append(sp.eye(X.shape[0]).tocsr())
    T_k.append(X)

    def chebyshev_recurrence(T_k_minus_one, T_k_minus_two, X):
        X_ = sp.csr_matrix(X, copy=True)
        return 2 * X_.dot(T_k_minus_one) - T_k_minus_two

    for i in range(2, k + 1):
        T_k.append(chebyshev_recurrence(T_k[-1], T_k[-2], X))

    for i in range(len(T_k)):
        T_k[i] = T_k[i].todense()

    return T_k
