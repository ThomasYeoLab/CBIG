import numpy as np

from scipy.spatial.distance import pdist, squareform
import scipy.sparse as sps


def distcorr(X, Y):
    """ Compute the distance correlation function

    >>> a = [1,2,3,4,5]
    >>> b = np.array([1,2,9,4,4])
    >>> np.allclose(distcorr(a, b), 0.762676242417)
    True
    """
    def allsame(x):
        return list(x) == list(x[::-1])
    if allsame(X) or allsame(Y):
        raise Exception("All elements of one input are equal, cannot divide by zero") 
    X = np.atleast_1d(X)
    Y = np.atleast_1d(Y)
    if np.prod(X.shape) == len(X):
        X = X[:, None]
    if np.prod(Y.shape) == len(Y):
        Y = Y[:, None]
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)
    n = X.shape[0]
    if Y.shape[0] != X.shape[0]:
        raise ValueError('Number of samples must match')
    a = squareform(pdist(X))
    b = squareform(pdist(Y))
    A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
    B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()

    dcov2_xy = (A * B).sum()/float(n * n)
    dcov2_xx = (A * A).sum()/float(n * n)
    dcov2_yy = (B * B).sum()/float(n * n)
    dcor = np.sqrt(dcov2_xy)/np.sqrt(np.sqrt(dcov2_xx) * np.sqrt(dcov2_yy))
    return dcor


def partial_distcorr(X, Y, Z):
    """ Compute the partial distance correlation function R(X,Y;Z)

    Reference: http://arxiv.org/abs/1310.2926

    >>> a = [1,2,3,4,5]
    >>> b = np.array([1,2,9,4,4])
    >>> c = np.array([9,5,6,7,8])
    >>> np.allclose(partial_distcorr(a, b, c), 0.80829037686547578)
    True

    >>> np.random.seed(0)
    >>> X = np.random.randn(30, 4)
    >>> Y  = np.random.randn(30, 4)
    >>> Z  = np.random.randn(30, 4)
    >>> np.allclose(partial_distcorr(X, Y, Z), 0.00042590779008772691)
    True

    >>> X = X + Z
    >>> Y = Y + Z
    >>> np.allclose(partial_distcorr(X, Y, Z), 0.20195314726364821)
    True
    """
    X = np.atleast_1d(X)
    Y = np.atleast_1d(Y)
    Z = np.atleast_1d(Z)
    if np.prod(X.shape) == len(X):
        X = X[:, None]
    if np.prod(Y.shape) == len(Y):
        Y = Y[:, None]
    if np.prod(Z.shape) == len(Z):
        Z = Z[:, None]
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)
    Z = np.atleast_2d(Z)
    n = X.shape[0]
    if Y.shape[0] != X.shape[0] or Y.shape[0] != Z.shape[0]:
        raise ValueError('Number of samples must match')
    a = squareform(pdist(X))
    b = squareform(pdist(Y))
    c = squareform(pdist(Z))
    # Compute U-centered matrics
    A = a - a.sum(axis=0)[None, :] / float(n - 2) - \
        a.sum(axis=1)[:, None]/float(n - 2) + a.sum() / float((n - 1) * (n - 2))
    B = b - b.sum(axis=0)[None, :] / float(n - 2) - \
        b.sum(axis=1)[:, None]/float(n - 2) + b.sum() / float((n - 1) * (n - 2))
    C = c - c.sum(axis=0)[None, :] / float(n - 2) - \
        c.sum(axis=1)[:, None]/float(n - 2) + c.sum() / float((n - 1) * (n - 2))
    A.flat[::(n+1)] = 0
    B.flat[::(n+1)] = 0
    C.flat[::(n+1)] = 0

    AB = (A * B).sum() #/float(n * (n - 3))
    AC = (A * C).sum() #/float(n * (n - 3))
    BC = (B * C).sum() #/float(n * (n - 3))
    nA = np.sqrt((A * A).sum())
    nB = np.sqrt((B * B).sum())
    nC = np.sqrt((C * C).sum())

    Rxy = AB / (nA * nB)
    Rxz = AC / (nA * nC)
    Ryz = BC / (nB * nC)

    pdcor = 0
    if not (np.allclose(Rxz, 1) or np.allclose(Ryz, 1)):
        pdcor = (Rxy - Rxz * Ryz) / (np.sqrt(1 - Rxz * Rxz) * np.sqrt(1 - Ryz * Ryz))
    else:
        PA = A - (AC / (nC * nC)) * C
        PB = B - (BC / (nC * nC)) * C
        nPA = np.sqrt((PA * PA).sum())
        nPB = np.sqrt((PB * PB).sum())
        if not np.allclose(nPA * nPB, 0):
            pdcor = (PA * PB).sum() / (nPA * nPB)
    return pdcor


def compute_nearest_neighbor_graph(K, n_neighbors=50):
    idx = np.argsort(K, axis=1)
    col = idx[:, -n_neighbors:].flatten()
    row = (np.array(range(K.shape[0]))[:, None] * np.ones((1, n_neighbors))).flatten().astype(int)
    A1 = sps.csr_matrix((np.ones((len(row))), (row, col)), shape=K.shape)
    A1 = (A1 + A1.transpose()) > 0
    idx1 = A1.nonzero()
    K = sps.csr_matrix((K.flat[idx1[0]*A1.shape[1] + idx1[1]],
                        A1.indices, A1.indptr))
    return K


def compute_affinity(X, method='markov', eps=None, metric='euclidean'):
    """Compute the similarity or affinity matrix between the samples in X

    :param X: A set of samples with number of rows > 1
    :param method: 'markov' or 'cauchy' kernel (default: markov)
    :param eps: scaling factor for kernel
    :param metric: metric to compute pairwise distances
    :return: a similarity matrix

    >>> X = np.array([[1,2,3,4,5], [1,2,9,4,4]])
    >>> np.allclose(compute_affinity(X, eps=1e3), [[1., 0.96367614], [ 0.96367614, 1.]])
    True

    >>> X = np.array([[1,2,3,4,5], [1,2,9,4,4]])
    >>> np.allclose(compute_affinity(X, 'cauchy', eps=1e3), [[0.001,  0.00096432], [ 0.00096432, 0.001 ]])
    True
    """
    import numpy as np
    D = squareform(pdist(X, metric=metric))
    if eps is None:
        k = int(max(2, np.round(D.shape[0] * 0.01)))
        eps = 2 * np.median(np.sort(D, axis=0)[k+1, :])**2
    if method == 'markov':
        affinity_matrix = np.exp(-(D * D) / eps)
    elif method == 'cauchy':
        affinity_matrix = 1./(D * D + eps)
    else:
        raise ValueError("Unknown method: {}".format(method))
    return affinity_matrix
