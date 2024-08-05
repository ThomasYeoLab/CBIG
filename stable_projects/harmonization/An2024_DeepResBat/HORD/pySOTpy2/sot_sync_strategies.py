#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
from __future__ import print_function
import numpy as np
import math
import logging
from poap.strategy import BaseStrategy, RetryStrategy

from HORD.pySOTpy2.experimental_design import SymmetricLatinHypercube, \
    LatinHypercube
from HORD.pySOTpy2.sampling_methods import CandidateDYCORS
from HORD.pySOTpy2.rbf_surfaces import CubicRBFSurface
from HORD.pySOTpy2.rbf_interpolant import RBFInterpolant
from HORD.pySOTpy2.utils import from_unit_box, round_vars, check_opt_prob
from HORD.pySOTpy2.rs_capped import RSPenalty

# Get module-level logger
logger = logging.getLogger(__name__)


class SyncStrategyNoConstraints(BaseStrategy):
    """Parallel synchronous optimization strategy without non-bound
    constraints.

    This class implements the parallel synchronous SRBF strategy
    described by Regis and Shoemaker.  After the initial experimental
    design (which is embarrassingly parallel), the optimization
    proceeds in phases.  During each phase, we allow nsamples
    simultaneous function evaluations.  We insist that these
    evaluations run to completion -- if one fails for whatever reason,
    we will resubmit it.  Samples are drawn randomly from around the
    current best point, and are sorted according to a merit function
    based on distance to other sample points and predicted function
    values according to the response surface.  After several
    successive significant improvements, we increase the sampling
    radius; after several failures to improve the function value, we
    decrease the sampling radius.  We restart once the sampling radius
    decreases below a threshold.
    """

    def __init__(self,
                 worker_id,
                 data,
                 response_surface,
                 maxeval,
                 nsamples,
                 exp_design=None,
                 sampling_method=None,
                 extra=None):
        """Initialize the optimization strategy.

        :param worker_id: Start ID in a multistart setting
        :param data: Problem parameter data structure
        :param response_surface: Surrogate model object
        :param maxeval: Function evaluation budget
        :param nsamples: Number of simultaneous fevals allowed
        :param exp_design: Experimental design
        :param sampling_method: Sampling method for finding
            points to evaluate
        :param extra: Points to be added to the experimental design
        """

        self.worker_id = worker_id
        self.data = data
        self.fhat = response_surface
        if self.fhat is None:
            self.fhat = RBFInterpolant(surftype=CubicRBFSurface, maxp=maxeval)

        self.maxeval = maxeval
        self.nsamples = nsamples
        self.extra = extra

        # Default to generate sampling points using Symmetric Latin Hypercube
        self.design = exp_design
        if self.design is None:
            if self.data.dim > 50:
                self.design = LatinHypercube(data.dim, data.dim + 1)
            else:
                self.design = SymmetricLatinHypercube(data.dim,
                                                      2 * (data.dim + 1))

        self.xrange = np.asarray(data.xup - data.xlow)

        # algorithm parameters
        self.sigma_min = 0.005
        self.sigma_max = 0.2
        self.sigma_init = 0.2

        self.failtol = max(5, data.dim)
        self.succtol = 3

        self.numeval = 0
        self.status = 0
        self.sigma = 0
        self.resubmitter = RetryStrategy()
        self.xbest = None
        self.fbest = np.inf
        self.fbest_old = None

        # Set up search procedures and initialize
        self.sampling = sampling_method
        if self.sampling is None:
            self.sampling = CandidateDYCORS(data)

        self.check_input()

        # Start with first experimental design
        self.sample_initial()

    def check_input(self):
        self.check_common()
        assert not hasattr(self.data, "eval_ineq_constraints"), \
            "Objective function has constraints,\n" \
            "SyncStrategyNoConstraints can't handle constraints"
        assert not hasattr(self.data, "eval_eq_constraints"), \
            "Objective function has constraints,\n" \
            "SyncStrategyNoConstraints can't handle constraints"

    def check_common(self):
        # Check evaluation budget
        if self.extra is None:
            assert self.maxeval >= self.design.npts, \
                "Experimental design is larger than the evaluation budget"
        else:
            assert self.maxeval >= self.design.npts + self.extra.shape[0], \
                "Experimental design + extra points exceeds the " + \
                "evaluation budget"
        # Check dimensionality
        assert self.design.dim == self.data.dim, \
            "Experimental design and optimization problem have " + \
            "different dimensions"
        if self.extra is not None:
            assert self.data.dim == self.extra.shape[1], \
                "Extra point and optimization problem have " + \
                "different dimensions"
        # Check that the optimization problem makes sense
        check_opt_prob(self.data)

    def proj_fun(self, x):
        x = np.atleast_2d(x)
        return round_vars(self.data, x)

    def log_completion(self, record):
        """Record a completed evaluation to the log.

        :param record: Record of the function evaluation
        :param penalty: Penalty for the given point
        """
        xstr = np.array_str(record.params[0],
                            max_line_width=np.inf,
                            precision=5,
                            suppress_small=True)
        if record.feasible:
            logger.info("{} {:.3e} @ {}".format("True", record.value, xstr))
        else:
            logger.info("{} {:.3e} @ {}".format("False", record.value, xstr))

    def adjust_step(self):
        """Adjust the sampling radius sigma.

        After succtol successful steps, we cut the sampling radius;
        after failtol failed steps, we double the sampling radius.

        :ivar Fnew: Best function value in new step
        :ivar fbest: Previous best function evaluation
        """
        # Initialize if this is the first adaptive step
        if self.fbest_old is None:
            self.fbest_old = self.fbest
            return

        # Check if we succeeded at significant improvement
        if self.fbest < self.fbest_old - 1e-3 * math.fabs(self.fbest_old):
            self.status = max(1, self.status + 1)
        else:
            self.status = min(-1, self.status - 1)
        self.fbest_old = self.fbest

        # Check if step needs adjusting
        if self.status <= -self.failtol:
            self.status = 0
            self.sigma /= 2
            logger.info("Reducing sigma")
        if self.status >= self.succtol:
            self.status = 0
            self.sigma = min([2.0 * self.sigma, self.sigma_max])
            logger.info("Increasing sigma")

    def sample_initial(self):
        """Generate and queue an initial experimental design.
        """
        if self.numeval == 0:
            logger.info("=== Start ===")
        else:
            logger.info("=== Restart ===")
        self.fhat.reset()
        self.sigma = self.sigma_init
        self.status = 0
        self.xbest = None
        self.fbest_old = None
        self.fbest = np.inf
        self.fhat.reset()

        start_sample = self.design.generate_points()
        assert start_sample.shape[1] == self.data.dim, \
            "Dimension mismatch between problem and experimental design"
        start_sample = from_unit_box(start_sample, self.data)
        if self.extra is not None:
            start_sample = np.vstack((start_sample, self.extra))

        for j in range(min(start_sample.shape[0],
                           self.maxeval - self.numeval)):
            start_sample[j, :] = self.proj_fun(
                start_sample[j, :])  # Project onto feasible region
            proposal = self.propose_eval(np.copy(start_sample[j, :]))
            self.resubmitter.rput(proposal)

        self.sampling.init(np.copy(start_sample), self.fhat,
                           self.maxeval - self.numeval)

    def sample_adapt(self):
        """Generate and queue samples from the search strategy
        """
        self.adjust_step()
        nsamples = min(self.nsamples, self.maxeval - self.numeval)
        new_points = self.sampling.make_points(npts=nsamples,
                                               xbest=np.copy(self.xbest),
                                               sigma=self.sigma,
                                               proj_fun=self.proj_fun)
        for i in range(nsamples):
            proposal = self.propose_eval(np.copy(np.ravel(new_points[i, :])))
            self.resubmitter.rput(proposal)

    def start_batch(self):
        """Generate and queue a new batch of points
        """
        if self.sigma < self.sigma_min:
            self.sample_initial()
        else:
            self.sample_adapt()

    def propose_action(self):
        """Propose an action
        """
        if self.numeval == self.maxeval:
            return self.propose_terminate()
        elif self.resubmitter.num_eval_outstanding == 0:
            self.start_batch()
        return self.resubmitter.get()

    def on_complete(self, record):
        """Handle completed function evaluation.

        When a function evaluation is completed we need to ask the constraint
        handler if the function value should be modified which is the case for
        say a penalty method. We also need to print the information to the
        logfile, update the best value found so far and notify the GUI that
        an evaluation has completed.

        :param record: Evaluation record
        """

        self.numeval += 1
        record.worker_id = self.worker_id
        record.worker_numeval = self.numeval
        record.feasible = True
        self.log_completion(record)
        self.fhat.add_point(np.copy(record.params[0]), record.value)
        if record.value < self.fbest:
            self.xbest = np.copy(record.params[0])
            self.fbest = record.value


class SyncStrategyPenalty(SyncStrategyNoConstraints):
    """Parallel synchronous optimization strategy with non-bound constraints.

    This is an extension of SyncStrategyNoConstraints that also works with
    bound constraints. We currently only allow inequality constraints, since
    the candidate based methods don't work well with equality constraints.
    We also assume that the constraints are cheap to evaluate, i.e., so that
    it is easy to check if a given point is feasible. More strategies that
    can handle expensive constraints will be added.

    We use a penalty method in the sense that we try to minimize:

    .. math::
        f(x) + \\mu \\sum_j (\\max(0, g_j(x))^2

    where :math:`g_j(x) \\leq 0` are cheap inequality constraints. As a
    measure of promising function values we let all infeasible points have
    the value of the feasible candidate point with the worst function value,
    since large penalties makes it impossible to distinguish between feasible
    points.

    When it comes to the value of :math:`\\mu`, just choose a very large value.


    """

    def __init__(self,
                 worker_id,
                 data,
                 response_surface,
                 maxeval,
                 nsamples,
                 exp_design=None,
                 sampling_method=None,
                 extra=None,
                 penalty=1e6):
        """Initialize the optimization strategy.

        :param worker_id: Start ID in a multistart setting
        :param data: Problem parameter data structure
        :param response_surface: Surrogate model object
        :param maxeval: Function evaluation budget
        :param nsamples: Number of simultaneous fevals allowed
        :param exp_design: Experimental design
        :param search_procedure: Search procedure for finding
            points to evaluate
        :param extra: Points to be added to the experimental design
        :param penalty: Penalty for violating constraints
        """

        # Evals wrapper for penalty method
        def penalty_evals(fhat, xx):
            penalty = self.penalty_fun(xx).T
            vals = fhat.evals(xx)
            if xx.shape[0] > 1:
                ind = (np.where(penalty <= 0.0)[0]).T
                if ind.shape[0] > 1:
                    ind2 = (np.where(penalty > 0.0)[0]).T
                    ind3 = np.argmax(np.squeeze(vals[ind]))
                    vals[ind2] = vals[ind3]
                    return vals
            return vals + penalty

        # Derivs wrapper for penalty method
        def penalty_derivs(fhat, xx):
            x = np.atleast_2d(xx)
            constraints = np.array(self.data.eval_ineq_constraints(x))
            dconstraints = self.data.deriv_ineq_constraints(x)
            constraints[np.where(constraints < 0.0)] = 0.0
            return np.atleast_2d(fhat.deriv(xx)) + \
                2 * self.penalty * np.sum(
                    constraints * np.rollaxis(dconstraints, 2), axis=2).T

        SyncStrategyNoConstraints.__init__(
            self, worker_id, data,
            RSPenalty(response_surface, penalty_evals, penalty_derivs),
            maxeval, nsamples, exp_design, sampling_method, extra)
        self.penalty = penalty

    def check_input(self):
        self.check_common()
        assert hasattr(
            self.data, "eval_ineq_constraints"), \
            "Objective function has no inequality constraints"
        assert not hasattr(self.data, "eval_eq_constraints"), \
            "Objective function has equality constraints,\n" \
            "SyncStrategyPenalty can't handle equality constraints"

    def penalty_fun(self, xx):
        """Computes the penalty for constraint violation

        :param xx: Points to compute the penalty for
        :return: Penalty for constraint violations
        """

        vec = np.array(self.data.eval_ineq_constraints(xx))
        vec[np.where(vec < 0.0)] = 0.0
        vec **= 2
        return self.penalty * np.asmatrix(np.sum(vec, axis=1))

    def on_complete(self, record):
        """Handle completed function evaluation.

        When a function evaluation is completed we need to ask the constraint
        handler if the function value should be modified which is the case for
        say a penalty method. We also need to print the information to the
        logfile, update the best value found so far and notify the GUI that
        an evaluation has completed.

        :param record: Evaluation record
        """
        x = np.zeros((1, record.params[0].shape[0]))
        x[0, :] = np.copy(record.params[0])
        penalty = self.penalty_fun(x)[0, 0]
        if penalty > 0.0:
            record.feasible = False
        else:
            record.feasible = True
        self.log_completion(record)
        self.numeval += 1
        record.worker_id = self.worker_id
        record.worker_numeval = self.numeval
        self.fhat.add_point(np.copy(record.params[0]), record.value)
        # Check if the penalty function is a new best
        if record.value + penalty < self.fbest:
            self.xbest = np.copy(record.params[0])
            self.fbest = record.value + penalty


class SyncStrategyProjection(SyncStrategyNoConstraints):
    """Parallel synchronous optimization strategy with non-bound constraints.
    It uses a supplied method to project proposed points onto the feasible
    region in order to always evaluate feasible points which is useful in
    situations where it is easy to project onto the feasible region and where
    the objective function is nonsensical for infeasible points.

    This is an extension of SyncStrategyNoConstraints that also works with
    bound constraints.
    """

    def __init__(self,
                 worker_id,
                 data,
                 response_surface,
                 maxeval,
                 nsamples,
                 exp_design=None,
                 sampling_method=None,
                 extra=None,
                 proj_fun=None):
        """Initialize the optimization strategy.

        :param worker_id: Start ID in a multistart setting
        :param data: Problem parameter data structure
        :param response_surface: Surrogate model object
        :param maxeval: Function evaluation budget
        :param nsamples: Number of simultaneous fevals allowed
        :param exp_design: Experimental design
        :param sampling_method: Search procedure for finding
            points to evaluate
        :param extra: Points to be added to the experimental design
        :param proj_fun: Projection operator
        """

        self.projection = proj_fun

        SyncStrategyNoConstraints.__init__(self, worker_id, data,
                                           response_surface, maxeval, nsamples,
                                           exp_design, sampling_method, extra)

    def check_input(self):
        self.check_common()
        assert hasattr(self.data, "eval_ineq_constraints") or \
            hasattr(self.data, "eval_eq_constraints"), \
            "Objective function has no constraints"

    def proj_fun(self, x):
        x = np.atleast_2d(x)
        for i in range(x.shape[0]):
            x[i, :] = self.projection(x[i, :])
        return x
