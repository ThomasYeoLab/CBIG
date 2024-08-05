#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
from copy import copy, deepcopy
import math
import numpy.linalg as la
from HORD.pySOTpy2.pyds import MassFunction


class EnsembleSurrogate:
    """Compute and evaluate an ensemble of interpolants.

    Maintains a list of surrogates and decides how to weights them
    by using Dempster-Shafer theory to assign pignistic probabilities
    based on statistics computed using LOOCV.

    :ivar nump: Current number of points
    :ivar maxp: Initial maximum number of points (can grow)
    :ivar rhs: Right hand side for interpolation system
    :ivar x: Interpolation points
    :ivar fx: Values at interpolation points
    :ivar dim: Number of dimensions
    :ivar model_list: List of surrogate models
    :ivar weights: Weight for each surrogate model
    :ivar surrogate_list: List of internal surrogate models for LOOCV
    """

    def __init__(self, model_list, maxp=100):
        self.nump = 0
        self.maxp = maxp
        self.x = None  # pylint: disable=invalid-name
        self.fx = None
        self.dim = None
        assert len(model_list) >= 2, "I need at least two models"
        self.model_list = model_list
        self.M = len(model_list)
        for i in range(self.M):
            self.model_list[i].reset()  # Models must be empty
        self.weights = None
        self.surrogate_list = None

    def reset(self):
        """Reset the interpolation."""
        self.nump = 0
        self.x = None
        self.fx = None
        for i in range(len(self.model_list)):
            self.model_list[i].reset()
        self.surrogate_list = None
        self.weights = None

    def _alloc(self, dim):
        """Allocate storage for x, fx, surrogate_list

        :param dim: Number of dimensions
        """
        maxp = self.maxp
        self.dim = dim
        self.x = np.zeros((maxp, dim))
        self.fx = np.zeros((maxp, 1))
        self.surrogate_list = [[None for _ in range(maxp)]
                               for _ in range(self.M)]

    def _realloc(self, dim, extra=1):
        """Expand allocation to accommodate more points (if needed)

        :param dim: Number of dimensions
        :param extra: Number of additional points to accommodate
        """
        if self.nump == 0:
            self._alloc(dim)
        elif self.nump + extra > self.maxp - 1:
            oldmaxp = self.maxp
            self.maxp = max([self.maxp * 2, self.maxp + extra])
            self.x.resize((self.maxp, dim))
            self.fx.resize((self.maxp, 1))
            # Expand the surrogate lists
            for i in range(self.M):
                for _ in range(self.maxp - oldmaxp):
                    self.surrogate_list[i].append(None)

    def _prob_to_mass(self, prob):
        """Internal method for building a mass function from probabilities

        :param prob: List of probabilities

        :return: Mass function
        """
        dictlist = []
        for i in range(len(prob)):
            dictlist.append([str(i + 1), prob[i]])
        return MassFunction(dict(dictlist))

    def _mean_squared_error(self, x, y):
        return np.sum((x - y)**2) / len(x)

    def _mean_abs_err(self, x, y):
        return np.sum(np.abs(x - y)) / len(x)

    def compute_weights(self):
        """Compute mode weights

        Given n observations we use n surrogates built with n-1 of the points
        in order to predict the value at the removed point. Based on these n
        predictions we calculate three different statistics:
            - Correlation coefficient with true function values
            - Root mean square deviation
            - Mean absolute error

        Based on these three statistics we compute the model weights by
        applying Dempster-Shafer theory to first compute the pignistic
        probabilities, which are taken as model weights.

        :return: Model weights
        """
        # Do the leave-one-out experiments
        loocv = np.zeros((self.M, self.nump))
        for i in range(self.M):
            for j in range(self.nump):
                loocv[i, j] = self.surrogate_list[i][j].eval(self.x[j, :])

        # Compute the model characteristics
        corr_coeff = np.ones(self.M)
        for i in range(self.M):
            corr_coeff[i] = np.corrcoef(
                np.vstack((loocv[i, :], self.get_fx().flatten())))[0, 1]

        root_mean_sq_err = np.ones(self.M)
        for i in range(self.M):
            root_mean_sq_err[i] = 1.0 / math.sqrt(
                self._mean_squared_error(self.get_fx().flatten(), loocv[i, :]))

        mean_abs_err = np.ones(self.M)
        for i in range(self.M):
            mean_abs_err[i] = 1.0 / self._mean_abs_err(self.get_fx().flatten(),
                                                       loocv[i, :])

        # Make sure no correlations are negative
        corr_coeff[np.where(corr_coeff < 0.0)] = 0.0
        if np.max(corr_coeff) == 0.0:
            corr_coeff += 1.0

        # Normalize the test statistics
        corr_coeff /= np.sum(corr_coeff)
        root_mean_sq_err /= np.sum(root_mean_sq_err)
        mean_abs_err /= np.sum(mean_abs_err)

        # Create mass functions based on the model characteristics
        m1 = self._prob_to_mass(corr_coeff)
        m2 = self._prob_to_mass(root_mean_sq_err)
        m3 = self._prob_to_mass(mean_abs_err)

        # Compute pignistic probabilities from Dempster-Shafer theory
        pignistic = m1.combine_conjunctive([m2, m3]).to_dict()
        self.weights = np.ones(self.M)
        for i in range(self.M):
            self.weights[i] = pignistic.get(str(i + 1))

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

        This function also updates the list of LOOCV surrogate models
        by cleverly
        just adding one point to n of the models. The scheme in which
        new models
        are built is illustrated below:

        2           1           1,2

        2,3         1,3         1,2         1,2,3

        2,3,4       1,3,4       1,2,4       1,2,3       1,2,3,4

        2,3,4,5     1,3,4,5     1,2,4,5     1,2,3,5     1,2,3,4     1,2,3,4,5

        :param xx: Point to add
        :param fx: The function value of the point to add
        """
        dim = len(xx)
        self._realloc(dim)
        self.x[self.nump, :] = xx
        self.fx[self.nump, :] = fx
        self.nump += 1
        # Update the leave-one-out models
        if self.nump == 2:
            for i in range(self.M):
                #  Add the first three models
                x0 = copy(self.x[0, :])
                x1 = copy(self.x[1, :])
                self.surrogate_list[i][0] = deepcopy(self.model_list[i])
                self.surrogate_list[i][0].add_point(x1, self.fx[1])
                self.surrogate_list[i][1] = deepcopy(self.model_list[i])
                self.surrogate_list[i][1].add_point(x0, self.fx[0])
                self.surrogate_list[i][2] = deepcopy(self.surrogate_list[i][1])
                self.surrogate_list[i][2].add_point(x1, self.fx[1])
        elif self.nump > 2:
            for i in range(self.M):
                for j in range(self.nump - 1):
                    self.surrogate_list[i][j].add_point(xx, fx)
                self.surrogate_list[i][self.nump] = deepcopy(
                    self.surrogate_list[i][self.nump - 1])
                self.surrogate_list[i][self.nump].add_point(xx, fx)
                # Point to the model with all points
                self.model_list[i] = self.surrogate_list[i][self.nump]
        self.weights = None

    def eval(self, xx, d=None):
        """Evaluate the interpolant at the point xx

        :param xx: Point where to evaluate

        :return: Value of the MARS interpolant at x
        """
        if self.weights is None:
            self.compute_weights()

        val = 0
        for i in range(self.M):
            val += self.weights[i] * self.model_list[i].eval(xx, d)
        return val

    def evals(self, xx, d=None):
        """Evaluate the MARS interpolant at the points xx

        :param xx: Points where to evaluate

        :return: Values of the MARS interpolant at x
        """
        if self.weights is None:
            self.compute_weights()

        vals = np.zeros((xx.shape[0], 1))
        for i in range(self.M):
            vals += self.weights[i] * self.model_list[i].evals(xx, d)

        return vals

    def deriv(self, x, d=None):
        """Evaluate the derivative of the MARS interpolant at x

        :param x: Data point

        :return: Derivative of the MARS interpolant at x
        """
        if self.weights is None:
            self.compute_weights()

        val = 0.0
        for i in range(self.M):
            val += self.weights[i] * self.model_list[i].deriv(x, d)
        return val


if __name__ == "__main__":

    from pySOT import RBFInterpolant
    from pySOT import CubicRBFSurface, TPSSurface, LinearRBFSurface

    fhat1 = RBFInterpolant(CubicRBFSurface, 1e-8, 100)
    fhat2 = RBFInterpolant(TPSSurface, 1e-8, 100)
    fhat3 = RBFInterpolant(LinearRBFSurface, 1e-8, 100)

    models = [fhat1, fhat2, fhat3]
    fhat = EnsembleSurrogate(models, 10)

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

    xs = np.random.rand(120, 2)
    for ii in range(100):
        xx = xs[ii, :]
        fx = test_f(xx)
        fhat.add_point(xx, fx)
    fhx = fhat.evals(xs[:5, :])
    print("Weights: " + str(fhat.weights))
    for ii in range(5):
        fx = test_f(xs[ii, :])
        print("Err: %e" % (abs(fx - fhx[ii]) / abs(fx)))
    for ii in range(10):
        xx = xs[100 + ii, :]
        fx = test_f(xx)
        dfx = test_df(xx)
        fhx = fhat.eval(xx)
        dfhx = fhat.deriv(xx)
        print("Err (interp): %e : %e" %
              (abs(fx - fhx) / abs(fx), la.norm(dfx - dfhx) / la.norm(dfx)))
