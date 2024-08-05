#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import numpy as np
import scipy.linalg as la
import scipy.spatial.distance as dist
from HORD.pySOTpy2.array_manager import ArrayManager


class BaseRBFSystem(object):
    """Base class for RBF system

    The RBF system is the set of centers and basis functions
    (polynomial and RBF), and the associated linear system.
    It does not involve actual function values or coefficients
    for interpolants.  One RBF system object may be shared among
    multiple interpolants.

    Derived classes should overload:
    - fill_basis - evaluate basis at points
    - basis_dot - compute linear combination of basis values
    - dbasis_dot - compute linear combination of basis gradients

    Derived classes may also overload the factor and solve routines.

    Attributes:
        new_points: Number of points since the last update
    """

    def __init__(self, dim, mtail, maxpts):
        """Initialize the system

        Args:
            dim: Dimension of the system
            mtail: Dimension of the polynomial tail
            maxpts: Initial estimate of maximum number of centers
        """
        self._mtail = mtail
        self._xmgr = ArrayManager((0, dim), (maxpts, dim))
        self._Mmgr = ArrayManager((self.nsys, self.nsys))
        self.new_points = 0

    @property
    def x(self):
        "Array of points"
        return self._xmgr.view

    @property
    def M(self):
        "System matrix"
        return self._Mmgr.view

    @property
    def dim(self):
        "Dimension of the ambient space"
        return self._xmgr.shape[1]

    @property
    def npt(self):
        "Number of centers"
        return self._xmgr.shape[0]

    @property
    def mtail(self):
        "Dimension of the polynomial tail space"
        return self._mtail

    @property
    def nsys(self):
        "System size"
        return self.npt + self.mtail

    @property
    def dirty(self):
        "Return true if there are new points"
        return self.new_points > 0

    def add_centers(self, x):
        """Add centers to the RBF interpolant

        Args:
            x: One or more new center locations (one per row if multiple)
        """
        if len(x.shape) == 1:
            self._xmgr.append(np.array([x]))
            self.new_points += 1
        else:
            self._xmgr.append(x)
            self.new_points += x.shape[0]

    def dist(self, y):
        """Compute distances of all points to y

        Args:
            y: One or more points (one per row if multiple)
        """
        if len(y.shape) == 1:
            return dist.cdist(self.x, np.expand_dims(y, axis=0)).ravel()
        else:
            return dist.cdist(self.x, y)

    def basis(self, y, d=None):
        """Compute the vector of basis functions at y

        Args:
            y: Point at which to evaluate basis functions
            d: Distance from y to RBF (not required)
        """
        b = np.zeros((self.nsys, 1))
        self.fill_basis(b, np.expand_dims(y, axis=0), d=None)
        b.shape = (self.nsys, )
        return b

    def compliance(self, y):
        """Get the compliance of the surface at y

        Args:
            y: Point at which to compute surface compliance
        """
        self.update()
        b = self.basis(y)
        return -np.dot(b, self.solve(b))

    def update(self):
        """Update the system matrix if points have been added
        """
        if self.new_points == 0:
            return
        x, nnew, Mmgr = self.x, self.new_points, self._Mmgr
        Mmgr.resize((self.nsys, self.nsys))
        M = Mmgr.view
        self.fill_basis(M[:, -nnew:], x[-nnew:, :])
        M[-nnew:, :] = M[:, -nnew:].T
        self.refactor()
        self.new_points = 0

    def refactor(self):
        """Compute factorization
        """
        eta = min(
            1e-5, 1e-16 * np.sqrt(
                la.norm(self.M, 1) * la.norm(self.M, np.inf))) * np.eye(
                    self.M.shape[0])
        self.lupiv = la.lu_factor(self.M + eta)

    def solve(self, rhs):
        """Solve with the system matrix

        Args:
            rhs: Right hand side in system

        Returns:
            solution to M*x = rhs
        """
        self.update()
        return la.lu_solve(self.lupiv, rhs)

    def fill_basis(self, B, y, d=None):
        """Compute a matrix of basis functions at y

        Args:
            B: Matrix to fill with basis function values
            y: Points at which to evaluate basis functions
            d: Distance from y to centers (optional)
        """
        pass

    def basis_dot(self, v, y, d=None):
        """Form a linear combination of basis functions

        Args:
            v: Coefficient vector in linear combination
            y: Point at which to evaluate basis functions
            d: Distance from y to centers (optional)
        """
        pass

    def dbasis_dot(self, v, y, d=None):
        """Form a linear combination of the basis gradients at y

        Args:
            v: Coefficient vector in linear combination
            y: Point at which to evaluate basis function gradients
            d: Distance from y to centers (optional)
        """
        pass


class SimpleRBFSystem(BaseRBFSystem):
    """Base class for RBF systems with simple polynomial tails

    Assumes that we use a low-degree polynomial tail (at most linear)
    with a standard basis.

    Derived classes should overload:
    - phi - evaluate RBF function
    - dphi_div_r - evaluate phi'(r)/r

    Attributes:
       dpoly: Degree bound for polynomial tail (derived class variable)
    """

    def __init__(self, dim, maxpts):
        """Initialize the system

        Args:
            dim: Dimension of the system
            dpoly: One more than the max poly degree (e.g. 2 for linear)
            maxpts: Initial estimate of maximum number of centers
        """
        dpoly = self.dpoly
        if dpoly == 0:
            mtail = 0
        elif dpoly == 1:
            mtail = 1
        elif dpoly == 2:
            mtail = dim + 1
        super(SimpleRBFSystem, self).__init__(dim, mtail, maxpts)

    def fill_basis(self, B, y, d=None):
        """Compute a matrix of basis functions at y

        Args:
            B: Matrix to fill with basis function values
            y: Points at which to evaluate basis functions
            d: Distance from y to centers (optional)
        """
        m = self.mtail
        if self.dpoly > 0:
            B[0, :] = 1
        if self.dpoly > 1:
            B[1:m, :] = y.T
        B[m:, :] = d if d is not None else self.dist(y)
        B[m:, :] = self.phi(B[m:, :])
        return B

    def basis_dot(self, v, y, d=None):
        """Form a linear combination of basis functions

        Args:
            v: Coefficient vector in linear combination
            y: Point at which to evaluate basis functions
            d: Distance from y to centers (optional)
        """

        if d is None:
            d = self.dist(y)
        m = self.mtail
        sy = np.dot(v[m:], self.phi(d))
        if self.dpoly > 0:
            sy += v[0]
        if self.dpoly > 1:
            sy += np.dot(y, v[1:m])
        return sy

    def dbasis_dot(self, v, y, d=None):
        """Form a linear combination of basis function gradients

        Args:
            v: Coefficient vector in linear combination
            y: Point at which to evaluate basis function gradients
            d: Distance from y to centers (optional)
        """
        if d is None:
            d = self.dist(y)
        x, m = self.x, self.mtail
        dsy = np.dot(self.dphi_div_r(d) * v[m:], y - x)
        if self.dpoly > 1:
            dsy += v[1:m]
        return dsy


class RBFSurface(object):
    """RBF surface

    An RBF surface builds an interpolant on top of an RBF system object.

    Attributes:
        rbfs: RBF system object
        coeff: interpolant coefficient vector
        dirty: flag whether there have been changes since coeff computation
    """

    def __init__(self, rbfs):
        """Initialize the surface

        Args:
            rbfs: RBF system object
        """
        self._fxmgr = ArrayManager((rbfs.nsys, ))
        self.coeff = None
        self.rbfs = rbfs
        self.dirty = True

    @property
    def x(self):
        "Array of points"
        return self.rbfs.x

    @property
    def fx(self):
        "Array of function values"
        return self._fxmgr.view

    @property
    def dim(self):
        "Dimension of the ambient space"
        return self.rbfs.dim

    @property
    def npt(self):
        "Number of centers"
        return self.rbfs.npt

    def set_fx(self, fx):
        "Replace f with transformed function values"
        self.fx[self.rbfs.mtail:] = fx
        self.dirty = True

    def fit(self):
        "Compute the fit"
        self.coeff = self.rbfs.solve(self.fx)
        self.dirty = False

    def update(self):
        "Update the fit if stale"
        if self.dirty:
            self.fit()

    def add_values(self, fx):
        """Add values to the system

        Args:
            fx: vector of values to append
        """
        self._fxmgr.append(fx)
        self.dirty = True

    def add_centers(self, x):
        """Add a center to the system

        Args:
            x: center to add
        """
        self.rbfs.add_centers(x)

    def add_points(self, x, fx):
        """Add points to the system

        NB: This function is only really safe if the RBF system is not
        shared.  To work with multiple surfaces associated with the same
        RBF system, add centers to the system object and values to the
        surfaces separately.

        Args:
            x: points to add to the RBF system
            fx: vector of values to append
        """
        self.add_centers(x)
        self.add_values(fx)

    def eval(self, y, d=None):
        """Evaluate the surface at y

        Args:
            y: Point or points at which to evaluate

        Returns:
            value of the surface at y
        """
        self.update()
        return self.rbfs.basis_dot(self.coeff, y, d)

    def deriv(self, y, d=None):
        """Evaluate the derivative at y

        Args:
            y: Point or points at which to evaluate

        Returns:
            gradient of the surface at y
        """
        self.update()
        return self.rbfs.dbasis_dot(self.coeff, y, d)

    def seminorm(self):
        "Get the seminorm of the surface"
        self.update()
        c, M = self.coeff, self.rbfs.M
        return np.dot(c, np.dot(M, c))

    def dseminorm(self, y, tau):
        """Compute change in seminorm associated with adding center at y

        Gutmann (2001) suggests choosing sample points in a global
        optimization by finding where the smallest seminorm change
        could be applied to the surface in order to achieve a given
        target value.

        Args:
            y: Point at which to add a hypothetical new center
            tau: Value of s(y) if augmented with hypothetical center

        Returns:
            g: Change in squared seminorm between original and updated surface
            dg: Gradient of g with respect to y

        """
        self.update()
        rbfs = self.rbfs
        d = rbfs.dist(y)
        b = rbfs.basis(y, d)
        w = rbfs.solve(b)
        mu = -np.dot(w, b)
        sy = np.dot(self.coeff, b)
        eta = (sy - tau) / mu
        g = (sy - tau) * eta
        dg = (2 * eta) * rbfs.dbasis_dot(self.coeff + eta * w, y, d)
        return g, dg

    def diseminorm(self, y, tau):
        """Compute 1/g(y) where g(y) is computed by dseminorm

        The function h(y) = -1/g(y) is potentially nicer to work
        with than just g(y), as it doesn't have any poles.

        Args:
            y: Point at which to add a hypothetical new center
            tau: Value of s(y) if augmented with hypothetical center

        Returns:
            h: Value of h(y) = 1/g(y)
            dh: Gradient of h with respect to y

        """
        self.update()
        rbfs = self.rbfs
        d = rbfs.dist(y)
        b = rbfs.basis(y, d)
        w = rbfs.solve(b)
        mu = -np.dot(w, b)
        sy = np.dot(self.coeff, b)
        gap = sy - tau
        h = mu / gap**2
        dh = -2 / gap**2 * rbfs.dbasis_dot(mu / gap * self.coeff + w, y, d)
        return h, dh


class SimpleRBFSurface(RBFSurface):
    """Simple RBF surface with a private RBF system object.

    Class variables:
        RBFSystem: Name of the RBF system object class

    Attributes:
        rbfs: RBF system object
        coeff: interpolant coefficient vector
        dirty: flag whether there have been changes since coeff computation
    """

    def __init__(self, x=None, fx=None, dim=None, maxpts=100):
        """Initialize the surface.

        Args:
            x: Initial points
            fx: Function values at f
            dim: Dimension of the space (required if x not given)
            maxpts: Estimate of max points (default: 100)
        """
        if x is not None:
            assert dim is None or dim == x.shape[1], "Dimension mismatch"
            assert fx is not None, "Must supply function values"
            dim = x.shape[1]
        rbfs = self.RBFSystem(dim, maxpts)
        super(SimpleRBFSurface, self).__init__(rbfs)
        if x is not None:
            self.add_points(x, fx)
