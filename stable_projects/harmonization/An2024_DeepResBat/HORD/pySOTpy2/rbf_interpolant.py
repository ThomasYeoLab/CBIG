#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
import numpy.linalg as la
from HORD.pySOTpy2.rbf_surfaces import CubicRBFSurface


class RBFInterpolant(object):
    """Compute and evaluate RBF interpolant.

    Manages an expansion of the form

    .. math::
        f(x) = \\sum_j c_j \\phi(\\|x-x_j\\|) + \\sum_j \\lambda_j p_j(x)

    where the functions :math:`p_j(x)` are low-degree polynomials.
    The fitting equations are

    .. math::

        \\begin{bmatrix} \\eta I & P^T \\\\ P & \\Phi+\\eta I \\end{bmatrix}
        \\begin{bmatrix} \\lambda \\\\ c \\end{bmatrix} =
        \\begin{bmatrix} 0 \\\\ f \\end{bmatrix}

    where :math:`P_{ij} = p_j(x_i)` and
    :math:`\\Phi_{ij}=\\phi(\\|x_i-x_j\\|)`.
    The regularization parameter :math:`\\eta` allows us to avoid problems
    with potential poor conditioning of the system.

    :ivar phi: Kernel function
    :ivar P: Tail functions
    :ivar dphi: Derivative of kernel function
    :ivar dP: Gradient of tail functions
    :ivar ntail: Number of tail functions
    :ivar nump: Current number of points
    :ivar maxp: Initial maximum number of points (can grow)
    :ivar A: Interpolation system matrix
    :ivar rhs: Right hand side for interpolation system
    :ivar x: Interpolation points
    :ivar fx: Values at interpolation points
    :ivar c: Expansion coefficients
    :ivar dim: Number of dimensions
    :ivar ntail: Number of tail functions
    """

    def __init__(self, surftype=CubicRBFSurface, maxp=100):
        self.maxp = maxp
        self.surface = None
        self.surftype = surftype

    @property
    def nump(self):
        return 0 if self.surface is None else self.surface.npt

    @property
    def x(self):
        return self.surface.x

    @property
    def fx(self):
        return self.surface.fx

    def reset(self):
        """Re-set the interpolation."""
        self.surface = None

    def get_x(self):
        """Get the list of data points

        :return: List of data points
        """
        return self.x

    def get_fx(self):
        """Get the list of function values for the data points.

        :return: List of function values
        """
        return np.expand_dims(self.fx, axis=1)

    def add_point(self, xx, fx):
        """Add a new function evaluation

        :param xx: Point to add
        :param fx: The function value of the point to add
        """
        if self.surface is None:
            self.surface = self.surftype(dim=len(xx), maxpts=self.maxp)
        else:
            assert xx.shape == (self.surface.dim, )
        self.surface.add_points(np.array([xx]), np.array([fx]))

    def transform_fx(self, fx):
        """Replace f with transformed function values for the fitting

        :param fx: Transformed function value
        """
        assert self.surface.npt == fx.shape[0]
        self.surface.set_fx(fx.ravel())

    def eval(self, xx, d=None):
        """Evaluate the rbf interpolant at the point xx

        :param xx: Point where to evaluate
        :return: Value of the rbf interpolant at x
        """
        assert xx.shape == (self.surface.dim, )
        return self.surface.eval(xx, d)

    def evals(self, xx, d=None):
        """Evaluate the rbf interpolant at the points xx

        :param xx: Points where to evaluate
        :return: Values of the rbf interpolant at x
        """
        assert xx.shape[1] == self.surface.dim
        return np.expand_dims(self.surface.eval(xx, d), axis=1)

    def deriv(self, x, d=None):
        """Evaluate the derivative of the rbf interpolant at x

        :param x: Data point
        :return: Derivative of the rbf interpolant at x
        """
        assert x.shape == (self.surface.dim, )
        return self.surface.deriv(x, d)


# ====================================================================


def _main():
    """Main test routine"""

    def test_f(x):
        """Test function"""
        fx = x[1] * np.sin(x[0]) + x[0] * np.cos(x[1])
        return fx

    def test_df(x):
        """Derivative of test function"""
        dfx = np.array([
            x[1] * np.cos(x[0]) + np.cos(x[1]),
            np.sin(x[0]) - x[0] * np.sin(x[1])
        ])
        return dfx

    fhat = RBFInterpolant(CubicRBFSurface, 1e-8, 20)
    xs = np.random.rand(120, 2)
    for i in range(100):
        xx = xs[i, :]
        fx = test_f(xx)
        fhat.add_point(xx, fx)
    fhx = fhat.evals(xs[:5, :])
    for i in range(5):
        fx = test_f(xs[i, :])
        print("Err: %e" % (abs(fx - fhx[i]) / abs(fx)))
    for i in range(10):
        xx = xs[100 + i, :]
        fx = test_f(xx)
        dfx = test_df(xx)
        fhx = fhat.eval(xx)
        dfhx = fhat.deriv(xx)
        print("Err (interp): %e : %e" %
              (abs(fx - fhx) / abs(fx), la.norm(dfx - dfhx) / la.norm(dfx)))


if __name__ == "__main__":
    _main()
