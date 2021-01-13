"""Align embeddings using Procrustes method
"""

import numpy as np


def match_coords(C1, C2):
    idx = []
    idxlist = range(C2.shape[0])
    for pt1 in C1:
        idxremain = np.setdiff1d(idxlist, idx)
        idxmatch = np.argsort(np.sum((C2[idxremain, :] - pt1) * (C2[idxremain, :] - pt1), axis=1))[0]
        idx.append(idxremain[idxmatch])
    return idx


def get_weight_matrix(Acoord, Bcoord, idx):
    d = np.sqrt(np.sum((Acoord - Bcoord[idx, :]) * (Acoord - Bcoord[idx, :]), axis=1))
    epsilon = max(np.median(d), 0.0001)
    W = np.diag(np.exp( - (d * d)/epsilon))
    return W


def iterative_alignment(embeddings, n_iters=1):
    target = embeddings[0]
    realigned = [target]
    xfms = []
    # first pass
    for i, embedding in enumerate(embeddings[1:]):
        u, s, v = np.linalg.svd(target.T.dot(embedding), full_matrices=False)
        xfms.append(v.T.dot(u.T))
        realigned.append(embedding.dot(xfms[-1]))

    # get mean target
    # 1. random sampling (doesn't cover all points)
    #    - allows more dense sampling
    #    - keeps coordinates from the real anatomical space
    # 2. mean across subjects
    # multi-pass
    for i in range(1, n_iters):
        index = []
        target = np.mean(realigned, axis=0).squeeze()
        target = np.array(target)
        realigned = []
        xfms = []
        for i, embedding in enumerate(embeddings):
            u, s, v = np.linalg.svd(target.T.dot(embedding), full_matrices=False)
            xfms.append(v.T.dot(u.T))
            realigned.append(embedding.dot(xfms[-1]))
    return realigned, xfms


def iterative_alignment_with_coords(embeddings, coords=None, n_iters=1, n_samples=0.1, use_mean=False):
    target = embeddings[0]
    if coords is None:
        targetcoords = None
        dummycoords = np.random.rand(target.shape[0], 3)
        W = np.eye(target.shape[0])
        idx = range(target.shape[0])
    else:
        targetcoords = coords[0]
    realigned = [target]
    # first pass
    for i, embedding in enumerate(embeddings[1:]):
        if targetcoords is not None:
            idx = match_coords(targetcoords, coords[i + 1])
            W = get_weight_matrix(targetcoords, coords[i + 1], idx)
        u, s, v = np.linalg.svd(target.T.dot(W.dot(embedding[idx, :])))
        xfm = v.T.dot(u.T)
        realigned.append(embedding.dot(xfm))

    # get mean target
    # 1. random sampling (doesn't cover all points)
    #    - allows more dense sampling
    #    - keeps coordinates from the real anatomical space
    # 2. mean across subjects
    # multi-pass
    for i in range(n_iters):
        index = []
        target = []
        targetcoords = []
        basecoords = []
        if use_mean:
            target = np.mean(realigned, axis=0).squeeze()
            if coords is None:
                targetcoords = dummycoords
            else:
                targetcoords = np.mean(coords, axis=0)
        else:
            for i, embedding in enumerate(realigned):
                index.append(np.random.permutation(embedding.shape[0])[:int(n_samples*embedding.shape[0])])
                target.extend(embedding[index[-1]].tolist())
                if coords is None:
                    targetcoords.extend(dummycoords[index[-1]])
                    basecoords.append(dummycoords)
                else:
                    targetcoords.extend(coords[i][index[-1]])
                    basecoords.append(coords[i])
        target = np.array(target)
        targetcoords = np.array(targetcoords)
        for i, embedding in enumerate(realigned):
            if coords is None:
                basecoords.append(dummycoords)
            else:
                basecoords.append(coords[i])
        realigned = []
        xfms = []
        for i, embedding in enumerate(embeddings):
            idx = match_coords(targetcoords, basecoords[i])
            W = get_weight_matrix(targetcoords, basecoords[i], idx)
            u, s, v = np.linalg.svd(target.T.dot(W.dot(embedding[idx, :])))
            xfms.append(v.T.dot(u.T))
            realigned.append(embedding.dot(xfms[-1]))
    return realigned, xfms
