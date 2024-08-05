#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import numpy as np
try:
    import matplotlib.pyplot as plt
    plotting_on = True
except BaseException:
    plotting_on = False
    pass


def to_unit_box(x, data):
    return (np.copy(x) - data.xlow) / (data.xup - data.xlow)


def from_unit_box(x, data):
    return data.xlow + (data.xup - data.xlow) * np.copy(x)


def unit_rescale(xx):
    """Shift and rescale elements of a vector to the unit interval
    """
    xmax = np.amax(xx)
    xmin = np.amin(xx)
    if xmax == xmin:
        return np.ones(xx.shape)
    else:
        return (xx - xmin) / (xmax - xmin)


def round_vars(data, x):
    """Round integer variables to closest integer
    """
    if len(data.integer) > 0:
        # Round the original ranged integer variables
        x[:, data.integer] = np.round(x[:, data.integer])
        # Make sure we don't violate the bound constraints
        for i in data.integer:
            ind = np.where(x[:, i] < data.xlow[i])
            x[ind, i] += 1
            ind = np.where(x[:, i] > data.xup[i])
            x[ind, i] -= 1
    return x


def check_opt_prob(obj):
    """Routine for checking that an implementation of the optimization problem
    follows the standard. This method checks everything, but can't make
    sure that the objective function and constraint methods return values
    of the correct type since this would involve actually evaluating the
    objective function which isn't feasible when the evaluations are
    expensive. If some test fails, an exception is raised through assert.

    :param obj: Objective function
    """
    assert hasattr(obj, "dim"), \
        "Problem dimension required"
    assert hasattr(obj, "xlow"), \
        "Numpy array of lower bounds required"
    assert isinstance(obj.xlow, np.ndarray), \
        "Numpy array of lower bounds required"
    assert hasattr(obj, "xup"), \
        "Numpy array of upper bounds required"
    assert isinstance(obj.xup, np.ndarray), \
        "Numpy array of upper bounds required"
    assert hasattr(obj, "integer"), \
        "Integer variables must be specified"
    if len(obj.integer) > 0:
        assert isinstance(obj.integer, np.ndarray), \
            "Integer variables must be specified"
    else:
        assert isinstance(obj.integer, np.ndarray) or \
            isinstance(obj.integer, list), \
            "Integer variables must be specified"
    assert hasattr(obj, "continuous"), \
        "Continuous variables must be specified"
    if len(obj.continuous) > 0:
        assert isinstance(obj.continuous, np.ndarray), \
            "Continuous variables must be specified"
    else:
        assert isinstance(obj.continuous, np.ndarray) or \
            isinstance(obj.continuous, list), \
            "Continuous variables must be specified"

    # Check for logical errors
    assert isinstance(obj.dim, int) and obj.dim > 0, \
        "Problem dimension must be a positive integer."
    assert (len(obj.xlow) == obj.dim and len(obj.xup) == obj.dim), \
        "Incorrect size for xlow and xup"
    assert all(obj.xlow[i] < obj.xup[i] for i in range(obj.dim)), \
        "Lower bounds must be below upper bounds."
    if len(obj.integer) > 0:
        assert np.amax(obj.integer) < obj.dim and np.amin(obj.integer) >= 0, \
            "Integer variable index can't exceed " \
            "number of dimensions or be negative"
    if len(obj.continuous) > 0:
        assert np.amax(obj.continuous) < obj.dim and \
            np.amin(obj.continuous) >= 0, \
            "Continuous variable index can't exceed " \
            "number of dimensions or be negative"
    assert len(np.intersect1d(obj.continuous, obj.integer)) == 0, \
        "A variable can't be both an integer and continuous"
    assert len(obj.continuous) + len(obj.integer) == obj.dim, \
        "All variables must be either integer or continuous"


def progress_plot(controller, title='', interactive=False):
    if not plotting_on:
        print("Failed to import matplotlib.pyplot, aborting....")
        return

    # Extract function values from the controller, ignoring crashed evaluations
    fvals = np.array(
        [o.value for o in controller.fevals if o.value is not None])

    plt.figure()
    if interactive:
        plt.ion()
    plt.plot(np.arange(0, fvals.shape[0]), fvals, 'bo')  # Points
    plt.plot(np.arange(0, fvals.shape[0]),
             np.minimum.accumulate(fvals),
             'r-',
             linewidth=4.0)  # Best value found

    # Set limits
    ymin = np.min(fvals)
    ymax = np.max(fvals)
    plt.ylim(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin))

    plt.xlabel('Evaluations')
    plt.ylabel('Function Value')
    plt.title(title)
    plt.show()
