"""Generate a diffusion map embedding
"""

import numpy as np
has_sklearn = True
try:
    import sklearn
except ImportError:
    has_sklearn = False


def compute_diffusion_map(L, alpha=0.5, n_components=None, diffusion_time=0,
                          skip_checks=False, overwrite=False,
                          eigen_solver=None, return_result=False):
    """Compute the diffusion maps of a symmetric similarity matrix

        L : matrix N x N
           L is symmetric and L(x, y) >= 0

        alpha: float [0, 1]
            Setting alpha=1 and the diffusion operator approximates the
            Laplace-Beltrami operator. We then recover the Riemannian geometry
            of the data set regardless of the distribution of the points. To
            describe the long-term behavior of the point distribution of a
            system of stochastic differential equations, we can use alpha=0.5
            and the resulting Markov chain approximates the Fokker-Planck
            diffusion. With alpha=0, it reduces to the classical graph Laplacian
            normalization.

        n_components: int
            The number of diffusion map components to return. Due to the
            spectrum decay of the eigenvalues, only a few terms are necessary to
            achieve a given relative accuracy in the sum M^t.

        diffusion_time: float >= 0
            use the diffusion_time (t) step transition matrix M^t

            t not only serves as a time parameter, but also has the dual role of
            scale parameter. One of the main ideas of diffusion framework is
            that running the chain forward in time (taking larger and larger
            powers of M) reveals the geometric structure of X at larger and
            larger scales (the diffusion process).

            t = 0 empirically provides a reasonable balance from a clustering
            perspective. Specifically, the notion of a cluster in the data set
            is quantified as a region in which the probability of escaping this
            region is low (within a certain time t).

        skip_checks: bool
            Avoid expensive pre-checks on input data. The caller has to make
            sure that input data is valid or results will be undefined.

        overwrite: bool
            Optimize memory usage by re-using input matrix L as scratch space.

        References
        ----------

        [1] https://en.wikipedia.org/wiki/Diffusion_map
        [2] Coifman, R.R.; S. Lafon. (2006). "Diffusion maps". Applied and
        Computational Harmonic Analysis 21: 5-30. doi:10.1016/j.acha.2006.04.006
    """

    import numpy as np
    import scipy.sparse as sps

    use_sparse = False
    if sps.issparse(L):
        use_sparse = True

    if not skip_checks:
        if has_sklearn:
            from sklearn.manifold.spectral_embedding_ import _graph_is_connected
            if not _graph_is_connected(L):
                raise ValueError('Graph is disconnected')
        else:
            raise ImportError('Checks require scikit-learn, but not found')

    ndim = L.shape[0]
    if overwrite:
        L_alpha = L
    else:
        L_alpha = L.copy()

    if alpha > 0:
        # Step 2
        d = np.array(L_alpha.sum(axis=1)).flatten()
        d_alpha = np.power(d, -alpha)
        if use_sparse:
            L_alpha.data *= d_alpha[L_alpha.indices]
            L_alpha = sps.csr_matrix(L_alpha.transpose().toarray())
            L_alpha.data *= d_alpha[L_alpha.indices]
            L_alpha = sps.csr_matrix(L_alpha.transpose().toarray())
        else:
            L_alpha = d_alpha[:, np.newaxis] * L_alpha 
            L_alpha = L_alpha * d_alpha[np.newaxis, :]

    # Step 3
    d_alpha = np.power(np.array(L_alpha.sum(axis=1)).flatten(), -1)
    if use_sparse:
        L_alpha.data *= d_alpha[L_alpha.indices]
    else:
        L_alpha = d_alpha[:, np.newaxis] * L_alpha

    M = L_alpha

    from scipy.sparse.linalg import eigs, eigsh
    if eigen_solver is None:
        eigen_solver = eigs

    # Step 4
    func = eigen_solver
    if n_components is not None:
        lambdas, vectors = func(M, k=n_components + 1)
    else:
        lambdas, vectors = func(M, k=max(2, int(np.sqrt(ndim))))
    del M

    if func == eigsh:
        lambdas = lambdas[::-1]
        vectors = vectors[:, ::-1]
    else:
        lambdas = np.real(lambdas)
        vectors = np.real(vectors)
        lambda_idx = np.argsort(lambdas)[::-1]
        lambdas = lambdas[lambda_idx]
        vectors = vectors[:, lambda_idx]

    return _step_5(lambdas, vectors, ndim, n_components, diffusion_time,
                   return_result)


def _step_5(lambdas, vectors, ndim, n_components, diffusion_time, return_result):
    """
    This is a helper function for diffusion map computation.

    The lambdas have been sorted in decreasing order.
    The vectors are ordered according to lambdas.

    """
    psi = vectors/vectors[:, [0]]
    diffusion_times = diffusion_time
    if diffusion_time == 0:
        diffusion_times = np.exp(1. -  np.log(1 - lambdas[1:])/np.log(lambdas[1:]))
        lambdas = lambdas[1:] / (1 - lambdas[1:])
    else:
        lambdas = lambdas[1:] ** float(diffusion_time)
    lambda_ratio = lambdas/lambdas[0]
    threshold = max(0.05, lambda_ratio[-1])

    n_components_auto = np.amax(np.nonzero(lambda_ratio > threshold)[0])
    n_components_auto = min(n_components_auto, ndim)
    if n_components is None:
        n_components = n_components_auto
    embedding = psi[:, 1:(n_components + 1)] * lambdas[:n_components][None, :]

    if return_result:
        result = dict(lambdas=lambdas, vectors=vectors,
                      n_components=n_components, diffusion_time=diffusion_times,
                      n_components_auto=n_components_auto)
        return embedding, result
    else:
        return embedding


def compute_diffusion_map_psd(
        X, alpha=0.5, n_components=None, diffusion_time=0, return_result=False):
    """
    This variant requires L to be dense, positive semidefinite and entrywise
    positive with decomposition L = dot(X, X.T).

    """
    from scipy.sparse.linalg import svds

    # Redefine X such that L is normalized in a way that is analogous
    # to a generalization of the normalized Laplacian.
    d = X.dot(X.sum(axis=0)) ** (-alpha)
    X = X * d[:, np.newaxis]

    # Decompose M = D^-1 X X^T
    # This is like
    # M = D^-1/2 D^-1/2 X (D^-1/2 X).T D^1/2
    # Substituting U = D^-1/2 X we have
    # M = D^-1/2 U U.T D^1/2
    # which is a diagonal change of basis of U U.T
    # which itself can be decomposed using svd.
    d = np.sqrt(X.dot(X.sum(axis=0)))
    U = X / d[:, np.newaxis]

    if n_components is not None:
        u, s, vh = svds(U, k=n_components+1, return_singular_vectors=True)
    else:
        k = max(2, int(np.sqrt(ndim)))
        u, s, vh = svds(U, k=k, return_singular_vectors=True)

    # restore the basis and the arbitrary norm of 1
    u = u / d[:, np.newaxis]
    u = u / np.linalg.norm(u, axis=0, keepdims=True)
    lambdas = s*s
    vectors = u

    # sort the lambdas in decreasing order and reorder vectors accordingly
    lambda_idx = np.argsort(lambdas)[::-1]
    lambdas = lambdas[lambda_idx]
    vectors = vectors[:, lambda_idx]

    return _step_5(lambdas, vectors, X.shape[0], n_components, diffusion_time,
                   return_result)


if has_sklearn:
    from sklearn.base import BaseEstimator
    import scipy.sparse as sps
    from sklearn.neighbors import kneighbors_graph

    class DiffusionMapEmbedding(BaseEstimator):
        """Diffusion map embedding for non-linear dimensionality reduction.

        Forms an affinity matrix given by the specified function and
        applies spectral decomposition to the corresponding graph laplacian.
        The resulting transformation is given by the value of the
        eigenvectors for each data point.

        Note : Laplacian Eigenmaps is the actual algorithm implemented here.

        Read more in the :ref:`User Guide <spectral_embedding>`.

        Parameters
        ----------

        diffusion_time : float
            Determines the scaling of the eigenvalues of the Laplacian

        alpha : float, optional, default: 0.5
            Setting alpha=1 and the diffusion operator approximates the
            Laplace-Beltrami operator. We then recover the Riemannian geometry
            of the data set regardless of the distribution of the points. To
            describe the long-term behavior of the point distribution of a
            system of stochastic differential equations, we can use alpha=0.5
            and the resulting Markov chain approximates the Fokker-Planck
            diffusion. With alpha=0, it reduces to the classical graph Laplacian
            normalization.

        n_components : integer, default: 2
            The dimension of the projected subspace.

        eigen_solver : {None, 'eigs' or 'eigsh'}
            The eigenvalue decomposition strategy to use.

        random_state : int, RandomState instance or None, optional, default: None
            A pseudo random number generator used for the initialization of the
            lobpcg eigenvectors.  If int, random_state is the seed used by the
            random number generator; If RandomState instance, random_state is the
            random number generator; If None, the random number generator is the
            RandomState instance used by `np.random`. Used when ``solver`` ==
            'amg'.

        affinity : string or callable, default : "nearest_neighbors"
            How to construct the affinity matrix.
             - 'nearest_neighbors' : construct affinity matrix by knn graph
             - 'rbf' : construct affinity matrix by rbf kernel
             - 'markov': construct affinity matrix by Markov kernel
             - 'cauchy': construct affinity matrix by Cauchy kernel
             - 'precomputed' : interpret X as precomputed affinity matrix
             - callable : use passed in function as affinity
               the function takes in data matrix (n_samples, n_features)
               and return affinity matrix (n_samples, n_samples).

        gamma : float, optional
            Kernel coefficient for pairwise distance (rbf, markov, cauchy)

        metric : string, optional
            Metric for scipy pdist function used to compute pairwise distances
            for markov and cauchy kernels

        n_neighbors : int, default : max(n_samples/10 , 1)
            Number of nearest neighbors for nearest_neighbors graph building.

        use_variant : boolean, default : False
            Use a variant requires L to be dense, positive semidefinite and
            entrywise positive with decomposition L = dot(X, X.T).

        n_jobs : int, optional (default = 1)
            The number of parallel jobs to run.
            If ``-1``, then the number of jobs is set to the number of CPU cores.

        Attributes
        ----------

        embedding_ : array, shape = (n_samples, n_components)
            Spectral embedding of the training matrix.

        affinity_matrix_ : array, shape = (n_samples, n_samples)
            Affinity_matrix constructed from samples or precomputed.

        References
        ----------

        - Lafon, Stephane, and Ann B. Lee. "Diffusion maps and coarse-graining: A
          unified framework for dimensionality reduction, graph partitioning, and
          data set parameterization." Pattern Analysis and Machine Intelligence,
          IEEE Transactions on 28.9 (2006): 1393-1403.
          https://doi.org/10.1109/TPAMI.2006.184

        - Coifman, Ronald R., and Stephane Lafon. Diffusion maps. Applied and
          Computational Harmonic Analysis 21.1 (2006): 5-30.
          https://doi.org/10.1016/j.acha.2006.04.006

        """

        def __init__(self, diffusion_time=0, alpha=0.5, n_components=2,
                     affinity="nearest_neighbors", gamma=None,
                     metric='euclidean', random_state=None, eigen_solver=None,
                     n_neighbors=None, use_variant=False, n_jobs=1):
            self.diffusion_time = diffusion_time
            self.alpha = alpha
            self.n_components = n_components
            self.affinity = affinity
            self.gamma = gamma
            self.metric = metric
            self.random_state = random_state
            self.eigen_solver = eigen_solver
            self.n_neighbors = n_neighbors
            self.use_variant = use_variant
            self.n_jobs = n_jobs

        @property
        def _pairwise(self):
            return self.affinity == "precomputed"

        def _get_affinity_matrix(self, X, Y=None):
            """Calculate the affinity matrix from data
            Parameters
            ----------
            X : array-like, shape (n_samples, n_features)
                Training vector, where n_samples is the number of samples
                and n_features is the number of features.

                If affinity is "precomputed"
                X : array-like, shape (n_samples, n_samples),
                Interpret X as precomputed adjacency graph computed from
                samples.

            Returns
            -------
            affinity_matrix, shape (n_samples, n_samples)
            """
            if self.affinity == 'precomputed':
                self.affinity_matrix_ = X
                return self.affinity_matrix_
            if self.affinity == 'nearest_neighbors':
                if sps.issparse(X):
                    warnings.warn("Nearest neighbors affinity currently does "
                                  "not support sparse input, falling back to "
                                  "rbf affinity")
                    self.affinity = "rbf"
                else:
                    self.n_neighbors_ = (self.n_neighbors
                                         if self.n_neighbors is not None
                                         else max(int(X.shape[0] / 10), 1))
                    self.affinity_matrix_ = kneighbors_graph(X, self.n_neighbors_,
                                                             include_self=True,
                                                             n_jobs=self.n_jobs)
                    # currently only symmetric affinity_matrix supported
                    self.affinity_matrix_ = 0.5 * (self.affinity_matrix_ +
                                                   self.affinity_matrix_.T)
                    return self.affinity_matrix_
            if self.affinity == 'rbf':
                self.gamma_ = (self.gamma
                               if self.gamma is not None else 1.0 / X.shape[1])
                self.affinity_matrix_ = rbf_kernel(X, gamma=self.gamma_)
                return self.affinity_matrix_
            if self.affinity in ['markov', 'cauchy']:
                from .dist import compute_affinity
                self.affinity_matrix_ = compute_affinity(X,
                                                         method=self.affinity,
                                                         eps=self.gamma,
                                                         metric=self.metric)
                return self.affinity_matrix_
            self.affinity_matrix_ = self.affinity(X)
            return self.affinity_matrix_

        def fit(self, X, y=None):
            """Fit the model from data in X.

            Parameters
            ----------
            X : array-like, shape (n_samples, n_features)
                Training vector, where n_samples is the number of samples
                and n_features is the number of features.

                If affinity is "precomputed"
                X : array-like, shape (n_samples, n_samples),
                Interpret X as precomputed adjacency graph computed from
                samples.

            Returns
            -------
            self : object
                Returns the instance itself.
            """

            from sklearn.utils import check_array, check_random_state
            X = check_array(X, ensure_min_samples=2, estimator=self)

            random_state = check_random_state(self.random_state)
            if isinstance(self.affinity, (str,)):
                if self.affinity not in set(("nearest_neighbors", "rbf",
                                             "markov", "cauchy",
                                             "precomputed")):
                    raise ValueError(("%s is not a valid affinity. Expected "
                                      "'precomputed', 'rbf', 'nearest_neighbors' "
                                      "or a callable.") % self.affinity)
            elif not callable(self.affinity):
                raise ValueError(("'affinity' is expected to be an affinity "
                                  "name or a callable. Got: %s") % self.affinity)

            affinity_matrix = self._get_affinity_matrix(X)
            if self.use_variant:
                self.embedding_ = compute_diffusion_map_psd(affinity_matrix,
                                                            alpha=self.alpha,
                                                            n_components=self.n_components,
                                                            diffusion_time=self.diffusion_time)
            else:
                self.embedding_ = compute_diffusion_map(affinity_matrix,
                                                        alpha=self.alpha,
                                                        n_components=self.n_components,
                                                        diffusion_time=self.diffusion_time,
                                                        eigen_solver=self.eigen_solver)
            return self

        def fit_transform(self, X, y=None):
            """Fit the model from data in X and transform X.

            Parameters
            ----------
            X : array-like, shape (n_samples, n_features)
                Training vector, where n_samples is the number of samples
                and n_features is the number of features.

                If affinity is "precomputed"
                X : array-like, shape (n_samples, n_samples),
                Interpret X as precomputed adjacency graph computed from
                samples.

            Returns
            -------
            X_new : array-like, shape (n_samples, n_components)
            """
            self.fit(X)
            return self.embedding_
