# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Miscellaneous numerical / torch helpers used across the individual-level
DELSSOME pipeline.

This module bundles four small helper areas:

* **Numpy helpers ported from the CBIG MATLAB library** — Pearson
  correlation (:func:`CBIG_corr`) and a numerically-stable arctanh
  (:func:`CBIG_StableAtanh`) used by Fisher-z averaging of FC matrices.
* **Torch port of the same helpers** — :func:`CBIG_corr_torch` mirrors
  :func:`CBIG_corr` on torch tensors so gradients / GPU acceleration
  propagate through it. Used in validity checks on sampled ``wEI``
  profiles.
* **FIC parameterisation** — :func:`parameterize_myelin_rsfc` converts
  the 10 parameterisation coefficients (intercept + myelin + RSFC
  gradient per ``wEE``, ``wEI``, ``sigma`` plus the global ``G``)
  optimised by CMA-ES into the ``3N+1`` per-region FIC parameter
  vector (Wang et al. 2019; Zhang et al. 2024).
* **Torch device / dtype defaults** — :func:`get_device` and
  :func:`set_torch_default` are called at the start of every entry
  point so that all tensors created downstream live on the same device
  and have the same dtype.
"""

import numpy as np
import scipy.stats as stats
import torch

from DELSSOME_indiv.CBIG_constants import DEFAULT_DTYPE


# -----------------------------------------------------------------------------
# Numpy helpers ported from the CBIG MATLAB library
# -----------------------------------------------------------------------------


def CBIG_corr(s_series, t_series=None, need_pvalue=False):
    """
    Pearson correlation matrix between rows of ``s_series`` (and ``t_series``).

    Mirrors the convention of the CBIG MATLAB library: features are along
    axis 0 (vertical) and labels along axis 1.

    Args:
        s_series (array): ``[features, labels]``.
        t_series (array, optional): ``[features, t_labels]``. If None,
            computes the self-correlation of ``s_series``.
        need_pvalue (bool): also return the two-sided p-values.

    Returns:
        ``corr_matrix [labels, t_labels]`` and optionally ``p_matrix``.
    """
    s_series = s_series - np.mean(s_series, axis=0)
    s_series = s_series / np.linalg.norm(s_series, axis=0)
    if t_series is None:
        corr_matrix = np.matmul(np.transpose(s_series), s_series)
    else:
        t_series = t_series - np.mean(t_series, axis=0)
        t_series = t_series / np.linalg.norm(t_series, axis=0)
        corr_matrix = np.matmul(np.transpose(s_series), t_series)
    if need_pvalue:
        n = s_series.shape[0]
        dist = stats.beta(n / 2 - 1, n / 2 - 1, loc=-1, scale=2)
        p_matrix = np.zeros_like(corr_matrix)
        for i in range(p_matrix.shape[0]):
            for j in range(p_matrix.shape[1]):
                p_matrix[i, j] = 2 * dist.cdf(-abs(corr_matrix[i, j]))
        return corr_matrix, p_matrix
    return corr_matrix


def CBIG_StableAtanh(x, ensure_real=True):
    """
    Numerically stable arctanh — clips inputs to ``[-1+eps, 1-eps]`` so that
    ``arctanh(±1)`` does not yield ``±inf``. Used by :func:`fisher_average`.
    """
    x[x > (1 - np.finfo(float).eps)] = 1 - np.finfo(float).eps
    x[x < (-1 + np.finfo(float).eps)] = -1 + np.finfo(float).eps
    x = np.arctanh(x)
    if not ensure_real:
        x[np.logical_and(np.logical_not(np.isreal(x)),
                         np.real(x) > 0)] = np.arctanh(1 - np.finfo(float).eps)
        x[np.logical_and(np.logical_not(np.isreal(x)),
                         np.real(x) < 0)] = np.arctanh(-1 + np.finfo(float).eps)
    return x


# -----------------------------------------------------------------------------
# Torch port of the CBIG numerical helpers
# -----------------------------------------------------------------------------


def CBIG_corr_torch(s_series, t_series=None):
    """
    Pearson correlation matrix between rows of ``s_series`` (and ``t_series``).

    Mirrors :func:`CBIG_corr` but operates on torch tensors so that
    gradients (and GPU acceleration) propagate through it. Used in
    validity checks on sampled ``wEI`` profiles.
    """
    s_series = s_series - torch.mean(s_series, dim=0)
    s_series = s_series / torch.linalg.norm(s_series, dim=0)
    if t_series is None:
        return torch.matmul(torch.transpose(s_series, 0, 1), s_series)
    t_series = t_series - torch.mean(t_series, dim=0)
    t_series = t_series / torch.linalg.norm(t_series, dim=0)
    return torch.matmul(torch.transpose(s_series, 0, 1), t_series)


# -----------------------------------------------------------------------------
# FIC parameterisation helpers
# -----------------------------------------------------------------------------


def parameterize_myelin_rsfc(myelin, rsfc_gradient, param_coef):
    """
    Convert 10 parameterisation coefficients into the 3N+1 FIC parameters.

    In the FIC model used here, the per-region parameters ``wEE``, ``wEI``
    and ``sigma`` are not optimised independently. Instead they are
    expressed as a linear combination of an intercept, the regional
    myelin map and the regional RSFC gradient (Wang et al. 2019;
    Zhang et al. 2024). CMA-ES then optimises the 10 coefficients of
    that parameterisation plus the global coupling ``G``.

    Args:
        myelin (Tensor): ``[N, 1]`` regional myelin map.
        rsfc_gradient (Tensor): ``[N, 1]`` regional RSFC gradient.
        param_coef (Tensor): ``[10, param_sets]`` parameterisation coefficients.

    Returns:
        Tensor: ``[3N+1, param_sets]`` per-region FIC parameters
        ``[wEE; wEI; G; sigma]``.
    """
    w_EE = (param_coef[0]
            + param_coef[1] * myelin
            + param_coef[2] * rsfc_gradient)
    w_EI = (param_coef[3]
            + param_coef[4] * myelin
            + param_coef[5] * rsfc_gradient)
    G = param_coef[6]
    sigma = (param_coef[7]
             + param_coef[8] * myelin
             + param_coef[9] * rsfc_gradient)
    return torch.vstack((w_EE, w_EI, G, sigma))


# -----------------------------------------------------------------------------
# Torch device / dtype defaults
# -----------------------------------------------------------------------------


def get_device():
    """Return ``torch.device('cuda')`` if a GPU is visible, else ``cpu``."""
    if torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")
    return device


def set_torch_default(device=None, dtype=None):
    """Set the default torch device + dtype globally and return them."""
    if device is None:
        device = get_device()
    torch.set_default_device(device)

    if dtype is None:
        dtype = DEFAULT_DTYPE
    torch.set_default_dtype(dtype)
    return device, dtype
