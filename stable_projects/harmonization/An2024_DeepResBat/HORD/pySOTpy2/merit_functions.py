#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''
import scipy.spatial as scp
import numpy as np
from HORD.pySOTpy2.utils import unit_rescale


def candidate_merit_weighted_distance(cand, npts=1):
    """Weighted distance merit function for the candidate points based methods
    :param cand: Candidate point object
    :param npts: Number of points selected for evaluation

    :return: Points selected for evaluation
    """

    new_points = np.ones((npts, cand.data.dim))

    for i in range(npts):
        ii = cand.next_weight
        weight = cand.weights[(ii + len(cand.weights)) % len(cand.weights)]
        merit = weight * cand.fhvals + \
            (1 - weight) * (1.0 - unit_rescale(cand.dmerit))

        merit[cand.dmerit < cand.dtol] = np.inf
        jj = np.argmin(merit)
        cand.fhvals[jj] = np.inf
        new_points[i, :] = cand.xcand[jj, :]

        # Update distances and weights
        ds = scp.distance.cdist(cand.xcand, np.atleast_2d(new_points[i, :]))
        cand.dmerit = np.minimum(cand.dmerit, ds)
        cand.next_weight += 1

    return new_points
