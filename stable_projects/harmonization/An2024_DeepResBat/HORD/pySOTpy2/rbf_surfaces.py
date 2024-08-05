#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
import scipy.linalg as la
from HORD.pySOTpy2.rbf_base import SimpleRBFSystem
from HORD.pySOTpy2.rbf_base import SimpleRBFSurface


class LinearRBFSystem(SimpleRBFSystem):
    """Linear RBF system"""

    dpoly = 1
    avgdist = False

    def phi(self, r, epsilon=None):
        return r

    def dphi_div_r(self, r, epsilon=None):
        return (r * 0 + 1) / r


class CubicRBFSystem(SimpleRBFSystem):
    """Cubic RBF system"""

    dpoly = 2
    avgdist = False

    def phi(self, r, epsilon=None):
        return r**3

    def dphi_div_r(self, r, epsilon=None):
        return 3 * r


class QuinticRBFSystem(SimpleRBFSystem):
    """Quintic RBF system"""

    dpoly = 2
    avgdist = False

    def phi(self, r, epsilon=None):
        return r**5

    def dphi_div_r(self, r, epsilon=None):
        return 5 * (r**3)


class ThinPlateRBFSystem(SimpleRBFSystem):
    """Thin plate RBF system"""

    dpoly = 2
    avgdist = False
    tiny = np.finfo(np.double).tiny

    def phi(self, r, epsilon=None):
        return r * r * np.log(r + self.tiny)

    def dphi_div_r(self, r, epsilon=None):
        return 2 * np.log(r + self.tiny) + 1


"""
class MultiQuadraticRBFSystem(SimpleRBFSystem):
    " Multi-Quadratic RBF "

    dpoly = 2
    avgdist = True

    def phi(self, r, epsilon):
        return np.sqrt((1.0/epsilon * r) ** 2 + 1)

    def dphi_div_r(self, r, epsilon):
        return 1.0/((epsilon ** 2) * np.sqrt(1 + (1.0/epsilon * r) ** 2))


class InverseMultiQuadraticRBFSystem(SimpleRBFSystem):
    " Inverse Multi-Quadratic RBF "

    dpoly = 2
    avgdist = True

    def phi(self, r, epsilon):
        return 1.0/np.sqrt((1.0 / epsilon * r) ** 2 + 1)

    def dphi_div_r(self, r, epsilon):
        return - 1.0/((epsilon ** 2) * (np.sqrt(1 +
        (1.0/epsilon * r) ** 2)) ** 3)


class GaussianRBFSystem(SimpleRBFSystem):
    " Gaussian RBF "

    dpoly = 2
    avgdist = True

    def phi(self, r, epsilon):
        return np.exp(-(1.0 / epsilon * r) ** 2)

    def dphi_div_r(self, r, epsilon):
        return - (2.0/(epsilon**2)) * np.exp(-(1.0 / epsilon * r) ** 2)
"""


class LinearRBFSurface(SimpleRBFSurface):
    RBFSystem = LinearRBFSystem


class CubicRBFSurface(SimpleRBFSurface):
    RBFSystem = CubicRBFSystem


class QuinticRBFSurface(SimpleRBFSurface):
    RBFSystem = QuinticRBFSystem


class TPSSurface(SimpleRBFSurface):
    RBFSystem = ThinPlateRBFSystem


"""
class MultiQuadSurface(SimpleRBFSurface):
    RBFSystem = MultiQuadraticRBFSystem


class InvMultiQuadSurface(SimpleRBFSurface):
    RBFSystem = InverseMultiQuadraticRBFSystem


class GaussianSurface(SimpleRBFSurface):
    RBFSystem = GaussianRBFSystem
"""


def toyf(x):
    if len(x.shape) == 1:
        return np.cos(2 * x[0] + np.cos(3 * x[1]))
    else:
        return np.cos(2 * x[:, 0] + np.cos(3 * x[:, 1]))


def fd_grad(f, x, h=1e-6):
    df = np.zeros(x.shape)
    xp = x.copy()
    for j in range(x.shape[0]):
        xp[j] = x[j] + h
        fp = f(xp)
        xp[j] = x[j] - h
        fm = f(xp)
        xp[j] = x[j]
        df[j] = (fp - fm) / 2 / h
    return df


def test_evals(RBFSurface=CubicRBFSurface):
    npts = 5
    x = np.random.random((npts, 2))
    fx = toyf(x)
    s = RBFSurface(x, fx)
    sx = s.eval(x)
    relerr = la.norm(sx - fx) / la.norm(fx)
    assert relerr < 1e-8, "Surface inconsistency: {0:.1e}".format(relerr)


def test_surface(RBFSurface=CubicRBFSurface):
    npts = 5
    x = np.random.random((npts, 2))
    y = np.random.random((2, ))
    s = RBFSurface(x, toyf(x))
    dg = s.deriv(y)
    dg_fd = fd_grad(s.eval, y)
    relerr = la.norm(dg - dg_fd) / la.norm(dg)
    assert relerr < 1e-8, "Surface inconsistency: {0:.1e}".format(relerr)


def test_seminorm(RBFSurface=CubicRBFSurface):
    npts = 5
    y = np.random.random((2, ))
    x0 = np.random.random((npts, 2))
    x1 = np.zeros((npts + 1, 2))
    x1[:-1, :] = x0
    x1[-1, :] = y
    s0 = RBFSurface(x0, toyf(x0))
    s1 = RBFSurface(x1, toyf(x1))
    gref = s1.seminorm() - s0.seminorm()
    g = s0.dseminorm(y, s1.fx[-1])[0]
    relerr = np.abs(g - gref) / gref
    assert relerr < 1e-8, "Seminorm inconsistency: {0:.1e}".format(relerr)


def test_seminorm1(RBFSurface=CubicRBFSurface):
    npts = 5
    x = np.random.random((npts, 2))
    y = np.random.random((2, ))
    s = RBFSurface(x, toyf(x))
    g, dg = s.dseminorm(y, toyf(y))
    ns0 = s.seminorm()
    s.add_points(y, toyf(y))
    ns1 = s.seminorm()
    gref = ns1 - ns0
    relerr = np.abs(g - gref) / gref
    assert relerr < 1e-8, "Seminorm inconsistency: {0:.1e}".format(relerr)


def test_gutmann(RBFSurface=CubicRBFSurface):
    npts = 5
    tau = -1
    x = np.random.random((npts, 2))
    y = np.random.random((2, ))
    s = RBFSurface(x, toyf(x))
    g, dg = s.dseminorm(y, tau)
    dg_fd = fd_grad(lambda y: s.dseminorm(y, tau)[0], y)
    relerr = la.norm(dg - dg_fd) / la.norm(dg)
    assert relerr < 1e-8, "Gutmann merit inconsistency: {0:.1e}".format(relerr)
    h, dh = s.diseminorm(y, tau)
    dh_fd = fd_grad(lambda y: s.diseminorm(y, tau)[0], y)
    relerr = la.norm(dh - dh_fd) / la.norm(dh)
    assert relerr < 1e-8, "Gutmann merit inconsistency: {0:.1e}".format(relerr)


if __name__ == "__main__":
    classes = [CubicRBFSurface, TPSSurface, LinearRBFSurface]
    for c in classes:
        test_evals(c)
        test_surface(c)
        test_seminorm(c)
        test_seminorm1(c)
        test_gutmann(c)
