#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
from pyKriging.krige import kriging


class KrigingInterpolant:
    """Compute and evaluate Kriging interpolant.

    :ivar nump: Current number of points
    :ivar maxp: Initial maximum number of points (can grow)
    :ivar x: Interpolation points
    :ivar fx: Function values at interpolation points
    :ivar k: Kriging model instance
    :ivar updated: Flag that indicates whether a Kriging
        model is up-to-date or not
    """

    def __init__(self, maxp=100):
        self.nump = 0
        self.maxp = maxp
        self.x = None
        self.fx = None
        self.dim = None
        self.k = None
        self.updated = False

    def reset(self):
        """Reset the Kriging interpolant
        """
        self.nump = 0
        self.x = None
        self.fx = None
        self.updated = False

    def _alloc(self, dim):
        """Allocate storage for x, fx and rhs.

        :param dim: Number of dimensions
        """
        maxp = self.maxp
        self.dim = dim
        self.x = np.zeros((maxp, dim))
        self.fx = np.zeros((maxp, 1))

    def _realloc(self, dim, extra=1):
        """Expand allocation to accommodate more points (if needed)

        :param dim: Number of dimensions
        :param extra: Number of additional points to accommodate
        """
        if self.nump == 0:
            self._alloc(dim)
        elif self.nump + extra > self.maxp:
            self.maxp = max(self.maxp * 2, self.maxp + extra)
            self.x.resize((self.maxp, dim))
            self.fx.resize((self.maxp, 1))

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
        dim = len(xx)
        self._realloc(dim)
        self.x[self.nump, :] = xx
        self.fx[self.nump, :] = fx
        self.nump += 1

        if self.k is None:
            self.k = kriging(self.x[:self.nump + 1, :],
                             self.fx[:self.nump + 1, :], self.maxp)
            self.k.train()

        # add point to kriging model
        self.k.addPoint(xx, fx)
        self.updated = False

    def eval(self, xx, d=None):
        """Evaluate the Kriging interpolant at the point xx

        :param xx: Point where to evaluate
        :return: Value of the Kriging interpolant at x
        """
        if self.updated is False:
            self.k.train()
        self.updated = True

        fx = self.k.predict(xx.ravel())
        return fx

    def evals(self, xx, d=None):
        """Evaluate the Kriging interpolant at the points xx

        :param xx: Points where to evaluate
        :return: Values of the Kriging interpolant at x
        """
        if self.updated is False:
            self.k.train()
        self.updated = True

        length = xx.shape[0]
        fx = np.zeros(shape=(length, 1))
        for i in range(length):
            fx[i, 0] = self.eval(np.asarray(xx[i]))
        return fx

    def deriv(self, x, d=None):
        """Evaluate the derivative of the rbf interpolant at x

        :param x: Data point
        :return: Derivative of the rbf interpolant at x
        """
        # FIXME, To be implemented
        raise NotImplementedError


# ====================================================================


def _main():
    """Main test routine"""

    def test_f(x):
        """Test function"""
        fx = x[1] * np.sin(x[0]) + x[0] * np.cos(x[1])
        return fx

    fhat = KrigingInterpolant(50)
    print("fhat.maxp: %i" % fhat.maxp)

    # Set up more points
    xs = np.random.rand(50, 2)
    for i in range(40):
        xx = xs[i, :]
        fx = test_f(xx)
        fhat.add_point(xx, fx)
    fhx = fhat.evals(xs[:10, :])
    for i in range(10):
        fx = test_f(xs[i, :])
        print("Err: %e" % (np.abs(fx - fhx[i]) / np.abs(fx)))


if __name__ == "__main__":
    _main()
