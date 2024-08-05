#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import math
import numpy as np
import scipy.spatial as scp
from scipy.optimize import minimize
import scipy.stats as stats
from HORD.pySOTpy2.merit_functions import candidate_merit_weighted_distance
from HORD.pySOTpy2.utils import unit_rescale
from HORD.pySOTpy2.heuristic_algorithms import GeneticAlgorithm as GA

# ========================= MultiSearch =======================


class MultiSampling(object):
    """ A collection of Sampling Methods and weights so that the user
        can use multiple sampling methods for the same optimization
        problem. This object keeps an internal list of proposed points
        in order to be able to compute the minimum distance from a point
        to all proposed evaluations. This list has to be reset each time
        the optimization algorithm restarts"""

    def __init__(self, strategy_list, cycle):
        if cycle is None:
            cycle = range(len(strategy_list))
        if (not all(isinstance(i, int) for i in cycle)) or \
                np.min(cycle) < 0 or \
                np.max(cycle) > len(strategy_list) - 1:
            raise ValueError("Incorrect cycle!!")
        self.sampling_strategies = strategy_list
        self.nstrats = len(strategy_list)
        self.cycle = cycle
        self.current_strat = 0
        self.proposed_points = None
        self.data = strategy_list[0].data
        self.fhat = None
        self.avoid = None
        self.budget = None
        self.n0 = None

    def init(self, start_sample, fhat, budget):
        """Initilize the sampling method by providing the points in the
        experimental design, the surrogate model, and the evaluation budget

        :param start_sample: Points in the experimental design
        :param fhat: Surrogate model
        :param budget: Evaluation budget
        """
        self.proposed_points = start_sample
        self.fhat = fhat
        self.n0 = start_sample.shape[0]
        for i in range(self.nstrats):
            self.sampling_strategies[i].init(self.proposed_points, fhat,
                                             budget)

    def remove_point(self, x):
        """Remove x from proposed_points. Useful if x was never evaluated.

        :param x: Point to be removed

        :return: True if points was removed, False otherwise
        """
        idx = np.sum(np.abs(self.proposed_points - x), axis=1).argmin()
        if np.sum(np.abs(self.proposed_points[idx, :] - x)) < 1e-10:
            self.proposed_points = np.delete(self.proposed_points, idx, axis=0)
            for i in range(self.nstrats):
                self.sampling_strategies[i].remove_point(x)
            return True
        return False

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        """Create new candidate points. This call is ignored by the
        optimization
        based search strategies.

        :param npts: Number of points to select
        :param xbest: Best solution found so far
        :param sigma: Current sampling radius w.r.t the unit box
        :param subset: Coordinates to perturb
        :param proj_fun: Routine for projecting infeasible points
        onto the feasible region
        :param merit: Merit function for selecting candidate points

        :return: Points selected for evaluation"""

        new_points = np.zeros((npts, self.data.dim))

        # Figure out what we need to generate
        npoints = np.zeros((self.nstrats, ), dtype=int)
        for i in range(npts):
            npoints[self.cycle[self.current_strat]] += 1
            self.current_strat = (self.current_strat + 1) % len(self.cycle)

        # Now generate the points from one strategy at the time
        count = 0
        for i in range(self.nstrats):
            if npoints[i] > 0:
                new_points[count:count + npoints[i], :] = \
                    self.sampling_strategies[i].make_points(
                        npts=npoints[i], xbest=xbest,
                        sigma=sigma, subset=subset,
                        proj_fun=proj_fun,
                        merit=merit)

                count += npoints[i]
                # Update list of proposed points
                for j in range(self.nstrats):
                    if j != i:
                        self.sampling_strategies[j].proposed_points = \
                            self.sampling_strategies[i].proposed_points

        return new_points


class CandidateSRBF(object):
    """This is an implementation of the candidate points method that is
    proposed in the first SRBF paper. Candidate points are generated
    by making normally distributed perturbations with stdDev sigma
    around the best solution

    :ivar data: Optimization object
    :ivar weights: Weights used in the merit function
    :ivar numcand: Number of candidate points to generate
    :ivar proposed_points: List of points proposed by any search strategy
        since the last restart
    """

    def __init__(self, data, numcand=None, weights=None):
        """Initialize the SRBF method

        :param data: Optimization object
        :param numcand: Number of candidate points to generate
        :param weights: Weights used for the merit function
        """
        self.data = data
        self.fhat = None
        self.xrange = self.data.xup - self.data.xlow
        self.dtol = 1e-3 * math.sqrt(data.dim)
        self.weights = weights
        if self.weights is None:
            self.weights = [0.3, 0.5, 0.8, 0.95]
        self.proposed_points = None
        self.dmerit = None
        self.xcand = None
        self.fhvals = None
        self.next_weight = 0
        self.numcand = numcand
        if self.numcand is None:
            self.numcand = min([5000, 100 * data.dim])
        self.budget = None
        self.n0 = None
        self.fhat = None

        # Check that the inputs make sense
        assert isinstance(self.numcand, int) and self.numcand > 0, \
            "The number of candidate points has to be a positive integer"
        assert (isinstance(
            self.weights, np.ndarray) or isinstance(
                self.weights, list)) and max(self.weights) <= 1 \
            and min(self.weights) >= 0, "Incorrect weights"

    def init(self, start_sample, fhat, budget):
        """Initialize the sampling method by providing the points in the
        experimental design, the surrogate model, and the evaluation budget

        :param start_sample: Points in the experimental design
        :param fhat: Surrogate model
        :param budget: Evaluation budget
        """
        self.proposed_points = start_sample
        self.n0 = start_sample.shape[0]
        self.budget = budget
        self.fhat = fhat

    def remove_point(self, x):
        """Remove x from the list of proposed points.
        Useful if x was never evaluated.

        :param x: Point to be removed

        :return: True if points was removed, False otherwise
        """
        idx = np.sum(np.abs(self.proposed_points - x), axis=1).argmin()
        if np.sum(np.abs(self.proposed_points[idx, :] - x)) < 1e-10:
            self.proposed_points = np.delete(self.proposed_points, idx, axis=0)
            return True
        return False

    def __generate_cand__(self, scalefactors, xbest, subset):
        self.xcand = np.ones((self.numcand, self.data.dim)) * xbest
        for i in subset:
            lower, upper = self.data.xlow[i], self.data.xup[i]
            ssigma = scalefactors[i]
            self.xcand[:, i] = stats.truncnorm.rvs((lower - xbest[i]) / ssigma,
                                                   (upper - xbest[i]) / ssigma,
                                                   loc=xbest[i],
                                                   scale=ssigma,
                                                   size=self.numcand)

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        """Create new candidate points based on the best solution
        and the current value of sigma.

        :param npts: Number of points to select
        :param xbest: Best solution found so far
        :param sigma: Current sampling radius w.r.t the unit box
        :param subset: Coordinates to perturb
        :param proj_fun: Routine for projecting infeasible points
        onto the feasible region
        :param merit: merit function for selecting candidate points

        :return: Points selected for evaluation
        """

        if subset is None:
            subset = np.arange(0, self.data.dim)
        scalefactors = sigma * self.xrange

        # Make sure that the scale factors are correct for
        # the integer variables (at least 1)
        ind = np.intersect1d(self.data.integer, subset)
        if len(ind) > 0:
            scalefactors[ind] = np.maximum(scalefactors[ind], 1.0)

        # Generate candidate points
        self.__generate_cand__(scalefactors, xbest, subset)
        if proj_fun is not None:
            self.xcand = proj_fun(self.xcand)

        dists = scp.distance.cdist(self.xcand, self.proposed_points)
        fhvals = self.fhat.evals(self.xcand)

        self.dmerit = np.amin(np.asmatrix(dists), axis=1)
        self.fhvals = unit_rescale(fhvals)

        xnew = merit(self, npts)
        self.proposed_points = np.vstack(
            (self.proposed_points, np.asmatrix(xnew)))
        return xnew


class CandidateUniform(CandidateSRBF):
    """Create Candidate points by sampling uniformly in the domain"""

    def __generate_cand__(self, scalefactors, xbest, subset):
        self.xcand = np.ones((self.numcand, self.data.dim)) * xbest
        self.xcand[:, subset] = np.random.uniform(self.data.xlow[subset],
                                                  self.data.xup[subset],
                                                  (self.numcand, len(subset)))


class CandidateDYCORS(CandidateSRBF):
    """This is an implementation of DyCORS method to generate
    candidate points. The DyCORS method only perturbs a subset
    of the dimensions when perturbing the best solution. The
    probability for a dimension to be perturbed decreases after
    each evaluation and is capped in order to guarantee
    global convergence."""

    def __init__(self, data, numcand=None, weights=None):
        """Initialize the DYCORS method

        :param data: Optimization object
        :param numcand:  Number of candidate points to generate"""

        CandidateSRBF.__init__(self, data, numcand=numcand, weights=weights)
        self.minprob = np.min([1.0, 1.0 / self.data.dim])
        assert data.dim > 1, "You can't use DYCORS on a 1d problem"

        def probfun(numevals, budget):
            if budget < 2:
                return 0
            return min([20.0 / data.dim, 1.0]) * \
                (1.0 - (np.log(numevals + 1.0) / np.log(budget)))

        self.probfun = probfun

    def __generate_cand__(self, scalefactors, xbest, subset):
        ddsprob = self.probfun(self.proposed_points.shape[0] - self.n0,
                               self.budget - self.n0)
        ddsprob = np.max([self.minprob, ddsprob])

        nlen = len(subset)
        ar = (np.random.rand(self.numcand, nlen) < ddsprob)
        ind = np.where(np.sum(ar, axis=1) == 0)[0]
        ar[ind, np.random.randint(0, nlen - 1, size=len(ind))] = 1

        self.xcand = np.ones((self.numcand, self.data.dim)) * xbest
        for i in range(nlen):
            lower, upper = self.data.xlow[i], self.data.xup[i]
            ssigma = scalefactors[subset[i]]
            ind = np.where(ar[:, i] == 1)[0]
            self.xcand[ind, subset[i]] = stats.truncnorm.rvs(
                (lower - xbest[subset[i]]) / ssigma,
                (upper - xbest[subset[i]]) / ssigma,
                loc=xbest[subset[i]],
                scale=ssigma,
                size=len(ind))


class CandidateDDS(CandidateDYCORS):
    """This is an implementation of DDS method to generate
    candidate points. Only a few candidate points are generated
    and the candidate point with the lowest value predicted
    by the surrogate model is selected. The DDS method only
    perturbs a subset of the dimensions when perturbing the
    best solution. The probability for a dimension to be
    perturbed decreases after each evaluation and is capped
    in order to guarantee global convergence."""

    def __init__(self, data, numcand=None, weights=None):
        """Initialize the DDS method

        :param data: Optimization object
        :param numcand:  Number of candidate points to generate"""

        CandidateDYCORS.__init__(self, data, numcand=numcand, weights=weights)
        self.weights = np.array([1.0])
        self.numcand = max([0.5 * data.dim, 2])

        def probfun(numevals, budget):
            return 1.0 - (np.log(numevals + 1.0) / np.log(budget))

        self.probfun = probfun

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        """Create new candidate points based on the best
        solution and the current value of sigma.

        :param npts: Number of points to select
        :param xbest: Best solution found so far
        :param sigma: Current sampling radius w.r.t the unit box
        :param subset: Coordinates to perturb
        :param proj_fun: Routine for projecting infeasible points
        onto the feasible region
        :param merit: merit function for selecting candidate points

        :return: Points selected for evaluation"""

        new_points = np.zeros((npts, self.data.dim))
        for i in range(npts):
            new_points[i, :] = CandidateDYCORS.make_points(self,
                                                           npts=1,
                                                           xbest=xbest,
                                                           sigma=0.2,
                                                           subset=subset,
                                                           proj_fun=proj_fun)
        return new_points


class CandidateSRBF_INT(CandidateSRBF):
    """Candidate points are generated by perturbing ONLY the discrete
    variables using the SRBF strategy
    """

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.integer) > 0:
            return CandidateSRBF.make_points(self,
                                             npts=npts,
                                             xbest=xbest,
                                             sigma=sigma,
                                             subset=self.data.integer,
                                             proj_fun=proj_fun,
                                             merit=merit)
        else:
            return CandidateSRBF.make_points(self,
                                             npts=npts,
                                             xbest=xbest,
                                             sigma=sigma,
                                             subset=self.data.continuous,
                                             proj_fun=proj_fun,
                                             merit=merit)


class CandidateDYCORS_INT(CandidateDYCORS):
    """Candidate points are generated by perturbing ONLY the discrete
    variables using the DYCORS strategy"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.integer) > 0:
            return CandidateDYCORS.make_points(self,
                                               npts=npts,
                                               xbest=xbest,
                                               sigma=sigma,
                                               subset=self.data.integer,
                                               proj_fun=proj_fun,
                                               merit=merit)
        else:
            return CandidateDYCORS.make_points(self,
                                               npts=npts,
                                               xbest=xbest,
                                               sigma=sigma,
                                               subset=self.data.continuous,
                                               proj_fun=proj_fun,
                                               merit=merit)


class CandidateDDS_INT(CandidateDDS):
    """Candidate points are generated by perturbing
    ONLY the discrete variables using the DDS strategy"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.integer) > 0:
            return CandidateDDS.make_points(self,
                                            npts=npts,
                                            xbest=xbest,
                                            sigma=sigma,
                                            subset=self.data.integer,
                                            proj_fun=proj_fun,
                                            merit=merit)
        else:
            return CandidateDDS.make_points(self,
                                            npts=npts,
                                            xbest=xbest,
                                            sigma=sigma,
                                            subset=self.data.continuous,
                                            proj_fun=proj_fun,
                                            merit=merit)


class CandidateUniform_INT(CandidateUniform):
    """Candidate points are generated by perturbing ONLY
    the discrete variables using the uniform perturbations"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.integer) > 0:
            return CandidateUniform.make_points(self,
                                                npts=npts,
                                                xbest=xbest,
                                                sigma=sigma,
                                                subset=self.data.integer,
                                                proj_fun=proj_fun,
                                                merit=merit)
        else:
            return CandidateUniform.make_points(self,
                                                npts=npts,
                                                xbest=xbest,
                                                sigma=sigma,
                                                subset=self.data.continuous,
                                                proj_fun=proj_fun,
                                                merit=merit)


class CandidateSRBF_CONT(CandidateSRBF):
    """Candidate points are generated by perturbing ONLY
    the continuous variables using the SRBF strategy"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.continuous) > 0:
            return CandidateSRBF.make_points(self,
                                             npts=npts,
                                             xbest=xbest,
                                             sigma=sigma,
                                             subset=self.data.continuous,
                                             proj_fun=proj_fun,
                                             merit=merit)
        else:
            return CandidateSRBF.make_points(self,
                                             npts=npts,
                                             xbest=xbest,
                                             sigma=sigma,
                                             subset=self.data.integer,
                                             proj_fun=proj_fun,
                                             merit=merit)


class CandidateDYCORS_CONT(CandidateDYCORS):
    """Candidate points are generated by perturbing ONLY
    the continuous variables using the DYCORS strategy"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.continuous) > 0:
            return CandidateDYCORS.make_points(self,
                                               npts=npts,
                                               xbest=xbest,
                                               sigma=sigma,
                                               subset=self.data.continuous,
                                               proj_fun=proj_fun,
                                               merit=merit)
        else:
            return CandidateDYCORS.make_points(self,
                                               npts=npts,
                                               xbest=xbest,
                                               sigma=sigma,
                                               subset=self.data.integer,
                                               proj_fun=proj_fun,
                                               merit=merit)


class CandidateDDS_CONT(CandidateDDS):
    """Candidate points are generated by perturbing
    ONLY the discrete variables using the DDS strategy"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.continuous) > 0:
            return CandidateDDS.make_points(self,
                                            npts=npts,
                                            xbest=xbest,
                                            sigma=sigma,
                                            subset=self.data.continuous,
                                            proj_fun=proj_fun,
                                            merit=merit)
        else:
            return CandidateDDS.make_points(self,
                                            npts=npts,
                                            xbest=xbest,
                                            sigma=sigma,
                                            subset=self.data.integer,
                                            proj_fun=proj_fun,
                                            merit=merit)


class CandidateUniform_CONT(CandidateUniform):
    """Candidate points are generated by perturbing ONLY
    the continuous variables using uniform perturbations"""

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=candidate_merit_weighted_distance):
        if len(self.data.continuous) > 0:
            return CandidateUniform.make_points(self,
                                                npts=npts,
                                                xbest=xbest,
                                                sigma=sigma,
                                                subset=self.data.continuous,
                                                proj_fun=proj_fun,
                                                merit=merit)
        else:
            return CandidateUniform.make_points(self,
                                                npts=npts,
                                                xbest=xbest,
                                                sigma=sigma,
                                                subset=self.data.integer,
                                                proj_fun=proj_fun,
                                                merit=merit)


# Optimization based strategies


class GeneticAlgorithm(object):

    def __init__(self, data):
        """Initialize the Genetic Algorithm method

        :param data: Optimization object
        """

        self.data = data
        self.fhat = None
        self.dtol = 1e-3 * math.sqrt(data.dim)
        self.proposed_points = None
        self.budget = None
        self.fhat = None

    def init(self, start_sample, fhat, budget):
        """Initialize the sampling method by providing the points in the
        experimental design, the surrogate model, and the evaluation budget

        :param start_sample: Points in the experimental design
        :param fhat: Surrogate model
        :param budget: Evaluation budget
        """

        self.proposed_points = start_sample
        self.budget = budget
        self.fhat = fhat

    def remove_point(self, x):
        """Remove x from the list of proposed points.
        Useful if x was never evaluated.

        :param x: Point to be removed

        :return: True if points was removed, False otherwise
        """

        idx = np.sum(np.abs(self.proposed_points - x), axis=1).argmin()
        if np.sum(np.abs(self.proposed_points[idx, :] - x)) < 1e-10:
            self.proposed_points = np.delete(self.proposed_points, idx, axis=0)
            return True
        return False

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=None):
        """Create new points to evaluate using the GA

        :param npts: Number of points to select
        :param xbest: Best solution found so far (ignored)
        :param sigma: Current sampling radius w.r.t the unit box (ignored)
        :param subset: Coordinates to perturb
        :param proj_fun: Routine for projecting infeasible points onto
        the feasible region
        :param merit: Merit function for selecting candidate points (ignored)

        :return: Points selected for evaluation"""

        new_points = np.zeros((npts, self.data.dim))
        for i in range(npts):
            ga = GA(self.fhat.evals,
                    self.data.dim,
                    self.data.xlow,
                    self.data.xup,
                    popsize=max([2 * self.data.dim, 100]),
                    ngen=100,
                    proj_fun=proj_fun)
            x_min, f_min = ga.optimize()

            dist = np.atleast_2d(
                np.min(scp.distance.cdist(self.proposed_points,
                                          np.atleast_2d(x_min)),
                       axis=1)).T

            x_new = x_min
            if np.min(dist) < self.dtol:
                # Perturb the best solution until we satisfy the tolerance
                d = 0.0
                x_new = None
                while d < self.dtol:
                    x_new = x_min + self.dtol * \
                        np.random.randn(1, self.data.dim)
                    x_new = np.maximum(x_new, self.data.xlow)
                    x_new = np.minimum(x_new, self.data.xup)
                    d = np.atleast_2d(
                        np.min(scp.distance.cdist(self.proposed_points,
                                                  np.atleast_2d(x_new)),
                               axis=1)).T.min()

            new_points[i, :] = x_new
            self.proposed_points = np.vstack(
                (self.proposed_points, np.asarray(x_new)))

        return new_points


class MultiStartGradient(object):
    """ A wrapper around the scipy.optimize implementations of box-constrained
        gradient based minimization.

    Attributes:
        usecand: Indicates that this method is NOT candidate based
        data: Optimization object
        fhat: Original response surface object.
        proposed_points: List of points proposed by any search strategy
                         since the last restart
        objfun: The merit function to minimize (needs to have a
        gradient available)
        bounds: n x 2 matrix with lower and upper bound constraints
        numrestarts: Number of random starting points
        method: What optimization method to use. The following
        options are available:
            - L-BFGS-B: Quasi-Newton method of
                        Broyden, Fletcher, Goldfarb, and Shanno (BFGS)
            - TNC:      Truncated Newton algorithm

        Note: SLSQP is supposed to work with bound constraints
        but for some reason it
              sometimes violates the constraints anyway.
    """

    def __init__(self, data, method='L-BFGS-B', num_restarts=30):
        """Initialize the Multi-Start Gradient object

        :param data: Optimization object
        :param method: Optimization method to use for the surrogate
        :param num_restarts: Number of restarts
        """

        self.data = data
        self.fhat = None
        self.bounds = np.zeros((self.data.dim, 2))
        self.bounds[:, 0] = self.data.xlow
        self.bounds[:, 1] = self.data.xup
        self.eval = None
        self.deriv = None
        self.dtol = 1e-3 * math.sqrt(data.dim)
        self.proposed_points = None
        self.budget = None
        self.num_restarts = num_restarts
        self.x_best = None
        if (method == 'TNC') or (method == 'L-BFGS-B'):
            self.method = method
        else:
            self.method = 'L-BFGS-B'

    def init(self, start_sample, fhat, budget):
        """Initialize the sampling method by providing the points in the
        experimental design, the surrogate model, and the evaluation budget

        :param start_sample: Points in the experimental design
        :param fhat: Surrogate model
        :param budget: Evaluation budget
        """

        self.proposed_points = start_sample
        self.budget = budget
        self.fhat = fhat

    def remove_point(self, x):
        """Remove x from the list of proposed points.
        Useful if x was never evaluated.

        :param x: Point to be removed

        :return: True if points was removed, False otherwise
        """

        idx = np.sum(np.abs(self.proposed_points - x), axis=1).argmin()
        if np.sum(np.abs(self.proposed_points[idx, :] - x)) < 1e-10:
            self.proposed_points = np.delete(self.proposed_points, idx, axis=0)
            return True
        return False

    def make_points(self,
                    npts,
                    xbest,
                    sigma,
                    subset=None,
                    proj_fun=None,
                    merit=None):
        """Create new points to evaluate using the Multi-Start Gradient method

        :param npts: Number of points to select
        :param xbest: Best solution found so far (ignored)
        :param sigma: Current sampling radius w.r.t the unit box (ignored)
        :param subset: Coordinates to perturb
        :param proj_fun: Routine for projecting infeasible points
        onto the feasible region
        :param merit: Merit for selecting candidate points (ignored)

        :return: Points selected for evaluation"""

        def eval(x):
            return self.fhat.eval(x).ravel()

        def deriv(x):
            return self.fhat.deriv(x).ravel()

        new_points = np.zeros((npts, self.data.dim))
        for j in range(npts):
            fvals = np.zeros(self.num_restarts)
            xvals = np.zeros((self.num_restarts, self.data.dim))
            dists = np.zeros(self.num_restarts)

            for i in range(self.num_restarts):
                if i == 0 and j == 0:
                    x0 = np.array(xbest)
                else:
                    x0 = np.random.uniform(self.data.xlow, self.data.xup)

                res = minimize(eval,
                               x0,
                               method=self.method,
                               jac=deriv,
                               bounds=self.bounds)

                # Compute the distance to the proposed points
                xx = np.atleast_2d(res.x)
                if proj_fun is not None:
                    xx = proj_fun(xx)
                dist = np.atleast_2d(
                    np.min(scp.distance.cdist(self.proposed_points, xx),
                           axis=1)).T

                fvals[i] = res.fun
                xvals[i, :] = xx
                dists[i] = dist.min()

            if dists.max() > self.dtol:
                x_new = None
                f_new = np.inf
                for i in range(self.num_restarts):
                    if dists[i] > self.dtol and fvals[i] < f_new:
                        x_new = xvals[i, :]
                        f_new = fvals[i]
            else:
                # Perturb the best point
                d = -1.0
                x_new = None
                while d < self.dtol:
                    x_new = xvals[
                        np.argmin(fvals), :] + self.dtol * np.random.randn(
                            1, self.data.dim)
                    x_new = np.maximum(x_new, self.data.xlow)
                    x_new = np.minimum(x_new, self.data.xup)
                    if proj_fun is not None:
                        x_new = proj_fun(x_new)
                    d = np.atleast_2d(
                        np.min(scp.distance.cdist(self.proposed_points,
                                                  np.atleast_2d(x_new)),
                               axis=1)).T.min()

            new_points[j, :] = x_new
            self.proposed_points = np.vstack(
                (self.proposed_points, np.asarray(x_new)))

        return new_points
