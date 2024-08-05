#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np


def extend_size(d, dmax):
    """Determine size to allocate on resize operation.

    Args:
        d: Desired length
        dmax: Current allocation size

    Returns:
        New allocation size
    """
    if d <= dmax:
        return dmax
    else:
        return 2 * d


class ArrayManager(object):
    """Manage extendible NumPy array storage.

    An ArrayManager handles copying and reallocation for an extensible
    storage block (as a NumPy array) attached to a NumPy view object.
    The manager does not ask as a proxy for the view.

    Attributes:
        view: View on the current active array
        shape: Current array shape
        maxdims: Maximum dimensions allowed without reallocation
    """

    def __init__(self, dims, maxdims=None, **kwargs):
        if maxdims is None:
            maxdims = dims
        self._init_kwargs = kwargs
        self._dims = dims
        self._data = np.zeros(maxdims, **self._init_kwargs)
        self._view = self._make_view()

    def _make_view(self):
        return self._data[tuple(slice(0, d) for d in self._dims)]

    @property
    def maxdims(self):
        return self._data.shape

    @property
    def shape(self):
        return self._dims

    @shape.setter
    def shape(self, newdims):
        self.resize(newdims)

    @property
    def view(self):
        return self._view

    def resize(self, dims):
        "Resize the array to the requested size"
        dnew = [
            extend_size(d, dmax) for d, dmax in zip(dims, self._data.shape)
        ]
        if dnew != self._data.shape:
            data = np.zeros(dnew, **self._init_kwargs)
            slices = tuple(slice(0, d) for d in self._dims)
            data[slices] = self._data[slices]
            self._data = data
            self._view = self._data[slices]
        self._dims = tuple(dims)
        self._view = self._make_view()

    def append(self, data, axis=0):
        "Add data along a given axis"
        if not isinstance(data, np.ndarray):
            data = np.array(data)
        if len(data.shape) == len(self._dims) - 1:
            data = np.expand_dims(data, axis=axis)
        daxis0 = self._dims[axis]
        daxis1 = daxis0 + data.shape[axis]
        newdims = list(self._dims)
        newdims[axis] = daxis1
        self.resize(newdims)
        slices = [slice(0, d) for d in self._dims]
        slices[axis] = slice(daxis0, daxis1)
        self._data[tuple(slices)] = data


def test_array_manager():
    "Test behavior of ArrayManager"

    mgr = ArrayManager((2, 3))
    v = mgr.view
    v[0, 0] = 1
    v[1, 2] = 2
    assert (v == np.array([[1, 0, 0], [0, 0, 2]])).all(), "Wrong initial view"
    mgr.shape = (3, 3)
    v2 = mgr.view
    assert (v2 == np.array([[1, 0, 0], [0, 0, 2],
                            [0, 0, 0]])).all(), "Wrong second view"

    assert mgr.maxdims[0] == 6, "Resize not as expected"
    assert mgr.maxdims[1] == 3, "Resize not as expected"

    mgr.append(np.array([1, 1, 1]))
    v3 = mgr.view
    assert (v3 == np.array([[1, 0, 0], [0, 0, 2], [0, 0, 0],
                            [1, 1, 1]])).all(), "Wrong third view"

    assert mgr.maxdims[0] == 6, "Resize not as expected"
    assert mgr.maxdims[1] == 3, "Resize not as expected"

    mgr.append(np.array([[1], [2], [3], [4]]), axis=1)
    v4 = mgr.view
    assert (v4 == np.array([[1, 0, 0, 1], [0, 0, 2, 2], [0, 0, 0, 3],
                            [1, 1, 1, 4]])).all(), "Wrong fourth view"

    assert mgr.maxdims[0] == 6, "Resize not as expected"
    assert mgr.maxdims[1] == 8, "Resize not as expected"

    mgr2 = ArrayManager((5, ), dtype=np.int)
    assert mgr2.view.dtype == np.int, "dtype transfer fails"

    mgr2.append(np.array([[1]]))
    assert mgr2.view.dtype == np.int, "dtype transfer fails"


if __name__ == "__main__":
    test_array_manager()
