import numpy as np
import scipy.stats as stats
from pySOT import CandidateSRBF


class BoundaryDYCORS(CandidateSRBF):
    def __init__(self, data, numcand=None, weights=None):
        """Initialize the DYCORS method

        :param data: Optimization object
        :param numcand:  Number of candidate points to generate"""

        CandidateSRBF.__init__(self, data, numcand=numcand, weights=weights)
        self.minprob = np.min([1.0, 1.0/self.data.dim])
        assert data.dim > 1, "You can't use DYCORS on a 1d problem"

        def probfun(numevals, budget):
            if budget < 2:
                return 0
            return min([20.0/data.dim, 1.0]) * (1.0 - (np.log(numevals + 1.0) / np.log(budget)))
        self.probfun = probfun

    def __generate_cand__(self, scalefactors, xbest, subset):
        ddsprob = self.probfun(self.proposed_points.shape[0] - self.n0, self.budget - self.n0)
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
            tmp = stats.norm.rvs(loc=xbest[subset[i]], scale=ssigma, size=len(ind))
            self.xcand[ind, subset[i]] = np.clip(tmp, lower, upper)
