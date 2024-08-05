#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
from HORD.pySOTpy2.utils import from_unit_box, to_unit_box


class RSCapped(object):
    """Cap adapter for RBF response surface.

    This adapter takes an existing response surface and replaces it
    with a modified version in which any function values above the
    median are replaced by the median value.

    :ivar model: Original response surface
    :ivar fvalues: Function values
    """

    def __init__(self, model, transformation=None):
        """Initialize the response surface adapter

        :param model: Original response surface object
        :param transformation: Function value transformation
        """
        self.needs_update = False
        self.transformation = transformation
        if self.transformation is None:

            def transformation(fvalues):
                medf = np.median(fvalues)
                fvalues[fvalues > medf] = medf
                return fvalues

            self.transformation = transformation
        self.model = model
        self.fvalues = np.zeros((model.maxp, 1))
        self.nump = 0
        self.maxp = model.maxp

    @property
    def x(self):
        return self.get_x()

    @property
    def fx(self):
        return self.get_fx()

    def reset(self):
        """Reset the capped response surface
        """
        self.model.reset()
        self.fvalues[:] = 0
        self.nump = 0

    def add_point(self, xx, fx):
        """Add a new function evaluation

        :param xx: Point to add
        :param fx: The function value of the point to add
        """
        if self.nump >= self.fvalues.shape[0]:
            self.fvalues.resize(2 * self.fvalues.shape[0], 1)
        self.fvalues[self.nump] = fx
        self.nump += 1
        self.needs_update = True
        self.model.add_point(xx, fx)

    def get_x(self):
        """Get the list of data points

        :return: List of data points
        """
        return self.model.get_x()

    def get_fx(self):
        """Get the list of function values for the data points.

        :return: List of function values
        """
        return self.model.get_fx()

    def eval(self, xx, d=None):
        """Evaluate the capped rbf interpolant at the point xx

        :param xx: Point where to evaluate
        :return: Value of the capped rbf interpolant at x
        """
        self._apply_transformation()
        return self.model.eval(xx, d)

    def evals(self, xx, d=None):
        """Evaluate the capped rbf interpolant at the points xx

        :param xx: Points where to evaluate
        :return: Values of the capped rbf interpolant at x
        """
        self._apply_transformation()
        return self.model.evals(xx, d)

    def deriv(self, xx, d=None):
        """Evaluate the derivative of the rbf interpolant at x

        :param x: Data point
        :return: Derivative of the rbf interpolant at x
        """
        self._apply_transformation()
        return self.model.deriv(xx, d)

    def _apply_transformation(self):
        """ Apply the cap to the function values.
        """
        fvalues = np.copy(self.fvalues[0:self.nump])
        self.model.transform_fx(self.transformation(fvalues))


class RSPenalty(object):
    """Cap adapter for RBF response surface.

    This adapter takes an existing response surface and replaces it
    with a modified version in which any function values above the
    median are replaced by the median value.

    :ivar model: Original response surface
    :ivar fvalues: Function values
    """

    def __init__(self, model, evals, derivs):
        """Initialize the response surface adapter

        :param model: Original response surface object
        """
        self.needs_update = False
        self.model = model
        self.fvalues = np.zeros((model.maxp, 1))
        self.nump = 0
        self.maxp = model.maxp
        self.eval_method = evals
        self.deriv_method = derivs

    @property
    def x(self):
        return self.get_x()

    @property
    def fx(self):
        return self.get_fx()

    def reset(self):
        """Reset the capped response surface
        """
        self.model.reset()
        self.fvalues[:] = 0
        self.nump = 0

    def add_point(self, xx, fx):
        """Add a new function evaluation

        :param xx: Point to add
        :param fx: The function value of the point to add
        """
        if self.nump >= self.fvalues.shape[0]:
            self.fvalues.resize(2 * self.fvalues.shape[0], 1)
        self.fvalues[self.nump] = fx
        self.nump += 1
        self.needs_update = True
        self.model.add_point(xx, fx)

    def get_x(self):
        """Get the list of data points

        :return: List of data points
        """
        return self.model.get_x()

    def get_fx(self):
        """Get the list of function values for the data points.

        :return: List of function values
        """
        return self.eval_method(self.model, self.model.get_x())[0, 0]

    def eval(self, xx, d=None):
        """Evaluate the capped rbf interpolant at the point xx

        :param xx: Point where to evaluate
        :return: Value of the capped rbf interpolant at x
        """
        return self.eval_method(self.model, np.atleast_2d(xx)).ravel()

    def evals(self, xx, d=None):
        """Evaluate the capped rbf interpolant at the points xx

        :param xx: Points where to evaluate
        :return: Values of the capped rbf interpolant at x
        """

        return self.eval_method(self.model, xx)

    def deriv(self, xx, d=None):
        """Evaluate the derivative of the rbf interpolant at x

        :param x: Data point
        :return: Derivative of the rbf interpolant at x
        """

        return self.deriv_method(self.model, xx)


class RSUnitbox(object):
    """Cap adapter for RBF response surface.

    This adapter takes an existing response surface and replaces it
    with a modified version in which any function values above the
    median are replaced by the median value.

    :ivar model: Original response surface
    :ivar fvalues: Function values
    """

    def __init__(self, model, data):
        """Initialize the response surface adapter

        :param model: Original response surface object
        """
        self.needs_update = False
        self.model = model
        self.fvalues = np.zeros((model.maxp, 1))
        self.nump = 0
        self.maxp = model.maxp
        self.data = data

    @property
    def x(self):
        return self.get_x()

    @property
    def fx(self):
        return self.get_fx()

    def reset(self):
        """Reset the capped response surface
        """
        self.model.reset()
        self.fvalues[:] = 0
        self.nump = 0

    def add_point(self, xx, fx):
        """Add a new function evaluation

        :param xx: Point to add
        :param fx: The function value of the point to add
        """
        if self.nump >= self.fvalues.shape[0]:
            self.fvalues.resize(2 * self.fvalues.shape[0], 1)
        self.fvalues[self.nump] = fx
        self.nump += 1
        self.needs_update = True
        self.model.add_point(to_unit_box(xx, self.data), fx)

    def get_x(self):
        """Get the list of data points

        :return: List of data points
        """
        return from_unit_box(self.model.get_x(), self.data)

    def get_fx(self):
        """Get the list of function values for the data points.

        :return: List of function values
        """
        return self.model.get_fx()

    def eval(self, xx, d=None):
        """Evaluate the capped rbf interpolant at the point xx

        :param xx: Point where to evaluate
        :return: Value of the capped rbf interpolant at x
        """
        return self.model.eval(to_unit_box(xx, self.data))

    def evals(self, xx, d=None):
        """Evaluate the capped rbf interpolant at the points xx

        :param xx: Points where to evaluate
        :return: Values of the capped rbf interpolant at x
        """

        return self.model.evals(to_unit_box(xx, self.data))

    def deriv(self, xx, d=None):
        """Evaluate the derivative of the rbf interpolant at x

        :param x: Data point
        :return: Derivative of the rbf interpolant at x
        """

        return self.model.deriv(to_unit_box(xx, self.data))
