#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
import numpy.linalg as la


class PolyRegression(object):
    """Compute and evaluate a polynomial regression surface.

    :ivar nump: Current number of points
    :ivar maxp: Initial maximum number of points (can grow)
    :ivar x: Interpolation points
    :ivar fx: Function evaluations of interpolation points
    :ivar dim: Number of dimensions
    :ivar bounds: Upper and lower bounds, one row per dimension
    :ivar basisp: Multi-indices representing terms in a tensor poly basis
        Each row is a list of dim indices indicating a polynomial degree
        in the associated dimension.
    """

    def __init__(self, bounds, basisp, maxp=100):
        """Initialize the objects

        :param bounds: a (dims, 2) array of lower and upper bounds
        in each coord
        :param basisp: a (nbasis, dims) array, where the ith basis function is
            prod_j L_basisp(i,j)(x_j), L_k = the degree k Legendre polynomial
        """
        self.nump = 0
        self.maxp = maxp
        self.x = None  # pylint: disable=invalid-name
        self.fx = None
        self.bounds = bounds
        self.basisp = basisp
        self.updated = False

    def reset(self):
        """Reset the interpolation."""
        self.nump = 0
        self.x = None
        self.fx = None
        self.updated = False

    # FIXME
    def poised_for_interpolation(self, X):
        return True

    @property
    def dim(self):
        return self.bounds.shape[0]

    def _normalize(self, x):
        """Normalize points to the box [-1,1]^d
        """
        xx = np.copy(x)
        for k in range(x.shape[1]):
            length = self.bounds[k, 0]
            u = self.bounds[k, 1]
            w = u - length
            xx[:, k] = (x[:, k] - length) / w + (x[:, k] - u) / w
        return xx

    def _alloc(self):
        """Allocate storage for x and fx.
        """
        maxp = self.maxp
        self.x = np.zeros((maxp, self.dim))
        self.fx = np.zeros((maxp, 1))

    def _realloc(self, extra=1):
        """Expand allocation to accommodate more points (if needed)

        :param extra: Number of additional points to accommodate
        """
        if self.nump == 0:
            self._alloc()
        elif self.nump + extra > self.maxp:
            self.maxp = max(self.maxp * 2, self.maxp + extra)
            self.x.resize((self.maxp, self.dim))
            self.fx.resize((self.maxp, 1))

    def _plegendre(self, x):
        """Evaluate basis functions.

        :param x: Coordinates (one per row)
        :return: Basis functions for each coordinate with shape (npts,nbasis)
        """
        s = self.basisp
        Px = legendre(x, np.max(s))
        Ps = np.ones((x.shape[0], s.shape[0]))
        for i in range(s.shape[0]):
            for j in range(s.shape[1]):
                Ps[:, i] *= Px[:, j, s[i, j]]
        return Ps

    def _dplegendre(self, x):
        """Evaluate basis function gradients.

        :param x: Coordinates (one per row)
        :return: Gradients for each coordinate with shape (npts,dim,nbasis)
        """
        s = self.basisp
        Px, dPx = dlegendre(x, np.max(s))
        dPs = np.ones((x.shape[0], x.shape[1], s.shape[0]))
        for i in range(s.shape[0]):
            for j in range(s.shape[1]):
                for k in range(x.shape[1]):
                    if k == j:
                        dPs[:, k, i] *= dPx[:, j, s[i, j]]
                    else:
                        dPs[:, k, i] *= Px[:, j, s[i, j]]
        return dPs

    def _fit(self):
        """Compute a least squares fit."""
        A = self._plegendre(self._normalize(self.get_x()))
        self.beta = la.lstsq(A, self.get_fx())[0]

    def _predict(self, xx):
        """Evaluate on response surface.
        """
        return np.dot(self._plegendre(self._normalize(xx)), self.beta)

    def _predict_deriv(self, xx):
        """Predict derivative.
        """
        dfx = np.dot(self._dplegendre(self._normalize(xx)), self.beta)
        for j in range(xx.shape[1]):
            dfx[:, j] /= (self.bounds[j, 1] - self.bounds[j, 0]) / 2
        return dfx

    def get_x(self):
        """Get the list of data points

        :return: List of data points
        """
        return self.x[:self.nump, :]

    def get_fx(self):
        """Get the list of function values for the data points.

        :return: List of function values
        """
        return self.fx[:self.nump, :]

    def add_point(self, xx, fx):
        """Add a new function evaluation

        :param xx: Point to add
        :param fx: The function value of the point to add
        """
        self._realloc()
        self.x[self.nump, :] = xx
        self.fx[self.nump, :] = fx
        self.nump += 1
        self.updated = False

    def eval(self, xx, d=None):
        """Evaluate the regression surface at point xx

        :param xx: Point where to evaluate
        :return: Value of the polynomial at xx
        """
        if self.updated is False:
            self._fit()
        self.updated = True

        xx = np.expand_dims(xx, axis=0)
        fx = self._predict(xx)
        return fx[0]

    def evals(self, xx, d=None):
        """Evaluate the regression surface at points xx

        :param xx: Points where to evaluate
        :return: Values of the polynomial at xx
        """
        if self.updated is False:
            self._fit()
        self.updated = True

        return np.atleast_2d(self._predict(xx))

    def deriv(self, x):
        """Evaluate the derivative of the regression surface at x

        :param x: Data point
        :return: Derivative of the polynomial at x
        """

        if self.updated is False:
            self._fit()
        self.updated = True

        x = np.expand_dims(x, axis=0)
        dfx = self._predict_deriv(x)
        return dfx[0]


def legendre(x, d):
    """Evaluate Legendre polynomials at all coordinates in x.

    :param x: Array of coordinates
    :param d: Max degree of polynomials
    :return: A x.shape-by-d array of Legendre polynomial values
    """
    x = np.array(x)
    s = x.shape + (d + 1, )
    x = np.ravel(x)
    P = np.zeros((x.shape[0], d + 1))
    P[:, 0] = 1
    if d > 0:
        P[:, 1] = x
    for n in range(1, d):
        P[:, n + 1] = ((2 * n + 1) * (x * P[:, n]) - n * P[:, n - 1]) / (n + 1)
    return P.reshape(s)


def dlegendre(x, d):
    """Evaluate Legendre polynomial derivatives at all coordinates in x.

    :param x: Array of coordinates
    :param d: Max degree of polynomials
    :return: x.shape-by-d arrays of Legendre polynomial values and derivatives
    """
    x = np.array(x)
    s = x.shape + (d + 1, )
    x = np.ravel(x)
    P = np.zeros((x.shape[0], d + 1))
    dP = np.zeros((x.shape[0], d + 1))
    P[:, 0] = 1
    if d > 0:
        P[:, 1] = x
        dP[:, 1] = 1
    for n in range(1, d):
        P[:, n + 1] = ((2 * n + 1) * (x * P[:, n]) - n * P[:, n - 1]) / (n + 1)
        dP[:, n + 1] = (
            (2 * n + 1) * (P[:, n] + x * dP[:, n]) - n * dP[:, n - 1]) / (
                n + 1)
    return P.reshape(s), dP.reshape(s)


def basis_base(n, testf):
    """Generate list of shape functions for a subset of a TP poly space.

    :param n: Dimension of the space
    :param testf: Return True if a given multi-index is in range
    :return: An N-by-n matrix with S(i,j) = degree of variable j in shape i
    """
    snext = np.zeros((n, ), dtype=np.int32)
    done = False

    # Follow carry chain through
    s = []
    while not done:
        s.append(snext.copy())
        done = True
        for i in range(n):
            snext[i] += 1
            if testf(snext):
                done = False
                break
            snext[i] = 0
    return np.array(s)


def basis_TP(n, d):
    """Generate list of shape functions for TP poly space.

    :param n: Dimension of the space
    :param d: Degree bound
    :return: An N-by-n matrix with S(i,j) = degree of variable j in shape i
           There are N = n^d shapes.
    """
    return basis_base(n, lambda s: np.all(s <= d))


def basis_TD(n, d):
    """Generate list of shape functions for TP poly space.

    :param n: Dimension of the space
    :param d: Degree bound
    :return: An N-by-n matrix with S(i,j) = degree of variable j in shape i
    """
    return basis_base(n, lambda s: np.sum(s) <= d)


def basis_HC(n, d):
    """Generate list of shape functions for HC poly space.

    :param n: Dimension of the space
    :param d: Degree bound
    :return: An N-by-n matrix with S(i,j) = degree of variable j in shape i
    """
    return basis_base(n, lambda s: np.prod(s + 1) <= d + 1)


def basis_SM(n, d):
    """Generate list of shape functions for SM poly space.

    :param n: Dimension of the space
    :param d: Degree bound
    :return: An N-by-n matrix with S(i,j) = degree of variable j in shape i
    """

    def fSM(p):
        return p if p < 2 else np.ceil(np.log2(p))

    def fSMv(s):
        f = 0
        for j in range(s.shape[0]):
            f += fSM(s[j])
        return f

    return basis_base(n, lambda s: fSMv(s) <= fSM(d))


# ====================================================================


def test_legendre1():
    npt = 1001
    x = np.linspace(-1, 1, npt)
    w = np.ones((npt, 1)) / (npt - 1)
    w[0] /= 2
    w[-1] /= 2
    P = legendre(x, 4)
    M = np.dot(P.T, w * P)
    E = M - np.diag(1 / (2 * np.arange(5) + 1))
    relerr = la.norm(E) / la.norm(M)
    assert relerr < 1e-4, "Test Legendre orthonormality"


def test_legendre2():
    npt = 1000
    x = np.linspace(-1, 1, npt)
    h = 2.0 / (npt - 1)
    P, dP = dlegendre(x, 4)
    assert np.max(
        np.abs(dP[1:-2, :] - (P[2:-1, :] - P[0:-3, :]) / 2 / h
               )) < 1e-4, "Test Legendre derivs vs finite difference"


def test_poly():
    bounds = np.array([[0, 1], [1, 1.5]])
    basisp = basis_TD(2, 2)
    surf = PolyRegression(bounds, basisp)

    def ref_poly(xy):
        x = xy[:, 0]
        y = xy[:, 1]
        return x**2 + y**2 / 3.0 + x * y / 4.0 + x / 5.0 + y / 6.0 + 1.0

    def dref_poly(xy):
        x = xy[:, 0]
        y = xy[:, 1]
        return np.array(
            [2.0 * x + y / 4.0 + 1 / 5.0, 2.0 * y / 3.0 + x / 4.0 + 1 / 6.0])

    def add_point(x, y):
        xy = np.array([[x, y]])
        surf.add_point(xy, ref_poly(xy))

    add_point(0, 1)
    add_point(1, 1.5)
    add_point(0, 1.4)
    add_point(0.7, 1.1)
    add_point(0.5, 1.3)
    add_point(0.9, 1.2)
    add_point(0.1, 1.25)

    xtest = np.array([[0.123, 1.456]])
    sxtest = surf.evals(xtest)
    fxtest = ref_poly(xtest)
    relerr = np.abs((sxtest - fxtest) / fxtest)
    assert relerr[0] < 1e-12

    dsxtest = surf.deriv(xtest[0, :])
    dfxtest = dref_poly(xtest)
    relerr = np.abs((dsxtest - dfxtest) / dfxtest)
    assert np.max(relerr) < 1e-12


if __name__ == "__main__":
    test_legendre1()
    test_legendre2()
    test_poly()
