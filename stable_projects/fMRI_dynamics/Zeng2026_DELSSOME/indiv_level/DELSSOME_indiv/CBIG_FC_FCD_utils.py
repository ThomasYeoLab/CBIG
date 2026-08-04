# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Functional connectivity (FC) and FC-dynamics (FCD) utilities, plus the FC/FCD
cost components used by the DELSSOME framework.

Only the functions actually exercised by the individual-level pipeline are
included. ComBat / ComBatLS harmonisation is not part of this release, which
works directly with raw per-subject FC and FCD CDFs (harmonisation happens
inside the GAMLSS fit instead; see Section 2.7 of the paper).
"""

import torch
import numpy as np
from DELSSOME_indiv.CBIG_misc_utils import CBIG_StableAtanh, set_torch_default


# --------------------------------------------------------------------------- #
# Averaging FCs across runs
# --------------------------------------------------------------------------- #
def fisher_average(corr_mat):
    """Fisher-z-transformed average over the first axis of ``corr_mat``.

    Args:
        corr_mat (ndarray): ``[number, ...]`` correlation matrices.
    Returns:
        ndarray: average correlation matrix.
    """
    z_corr_mat = CBIG_StableAtanh(corr_mat)
    z_corr_ave = np.nanmean(z_corr_mat, axis=0)
    corr_ave = np.tanh(z_corr_ave)
    return corr_ave


def torch_corr_3D(vec_3d):
    """
    Pearson correlation between two stacked vectors per "M" sample.

    Args:
        vec_3d (Tensor): ``[M, 2, len]`` — the middle axis indexes (X, Y).
    Returns:
        Tensor: ``[M]`` correlations between X and Y for every M-th pair.
    """
    ex_ey = torch.mean(vec_3d, dim=-1)
    # ``unbiased=False`` matches ``torch.mean`` (denominator = N, not N - 1).
    std = torch.std(vec_3d, dim=-1, unbiased=False)
    xy = vec_3d[:, 0, :] * vec_3d[:, 1, :]
    e_xy = torch.mean(xy, dim=-1)
    cov = e_xy - ex_ey[:, 0] * ex_ey[:, 1]
    corr = cov / (std[:, 0] * std[:, 1])
    return corr


# --------------------------------------------------------------------------- #
# FC vectorisation helpers
# --------------------------------------------------------------------------- #
def create_fc_mask(N):
    """Boolean ``[N, N]`` mask for the strict upper triangle (off-diagonal)."""
    mask = torch.ones(N, N, dtype=torch.bool)
    mask = torch.triu(mask, 1)
    return mask


def vectorize_fc(fc_mat):
    """Extract the strict upper-triangle of an FC matrix (or batch thereof)."""
    if fc_mat.ndim == 2:
        N = fc_mat.shape[0]
        mask = create_fc_mask(N)
        return fc_mat[mask]
    if fc_mat.ndim == 3:
        _, N, _ = fc_mat.shape
        mask = create_fc_mask(N)
        return fc_mat[:, mask]
    raise ValueError(f"FC matrix must be 2D or 3D, got {fc_mat.ndim}D")


def devectorize_fc(fc_vec, N):
    """Inverse of :func:`vectorize_fc` — rebuild a symmetric FC matrix."""
    set_torch_default()
    if not isinstance(fc_vec, torch.Tensor):
        fc_vec = torch.tensor(fc_vec)

    expected_vec_len = N * (N - 1) // 2
    if fc_vec.ndim == 1:
        if fc_vec.shape[0] != expected_vec_len:
            raise ValueError(
                f"Vector length {fc_vec.shape[0]} doesn't match {expected_vec_len} for N={N}"
            )
        fc_mat = torch.zeros(N, N)
        mask = create_fc_mask(N)
        fc_mat[mask] = fc_vec
        fc_mat = fc_mat + fc_mat.T
        fc_mat.fill_diagonal_(1.0)
        return fc_mat

    if fc_vec.ndim == 2:
        M = fc_vec.shape[0]
        if fc_vec.shape[1] != expected_vec_len:
            raise ValueError(
                f"Vector length {fc_vec.shape[1]} doesn't match {expected_vec_len} for N={N}"
            )
        fc_mat = torch.zeros(M, N, N)
        mask = create_fc_mask(N)
        fc_mat[:, mask] = fc_vec
        fc_mat = fc_mat + fc_mat.transpose(-2, -1)
        for i in range(N):
            fc_mat[:, i, i] = 1.0
        return fc_mat
    raise ValueError(f"FC vector must be 1D or 2D, got {fc_vec.ndim}D")


# --------------------------------------------------------------------------- #
# FC / FCD computation from BOLD
# --------------------------------------------------------------------------- #
def calc_FC(bold: torch.Tensor):
    """Compute the FC matrix for each of M sets of BOLD signals.

    Args:
        bold (Tensor): ``[N, M, t_len]`` or ``[N, t_len]``.
    Returns:
        Tensor: FC matrix ``[M, N, N]`` (or ``[N, N]`` if M was singleton).
    """
    if bold.ndim == 2:
        bold = bold.unsqueeze(1)

    N = bold.shape[0]
    M = bold.shape[1]
    fc_mat = torch.zeros((M, N, N))
    for i in range(M):
        fc_mat[i] = torch.corrcoef(bold[:, i, :])

    if fc_mat.shape[0] == 1:
        fc_mat = fc_mat.squeeze(0)
    return fc_mat


def calc_FCD(bold: torch.Tensor,
             window_size=83,
             bins=10000,
             get_FCD_matrix=False):
    """
    Sliding-window FCD computation (Hansen et al., 2015; Zalesky et al., 2014).

    For each parameter set, FC is computed inside every length-``window_size``
    window of the BOLD signal; the Pearson correlation between these
    sliding-window FCs yields an ``[M, window_num, window_num]`` FCD matrix.
    Its off-diagonal entries are then histogrammed over ``bins`` bins from
    ``-1`` to ``1`` to obtain a discrete FCD distribution that can be matched
    against an empirical FCD CDF via the Kolmogorov-Smirnov distance.

    Args:
        bold (Tensor): ``[N, M, t_len]`` or ``[N, t_len]``.
        window_size (int): sliding-window length in time-points (~60 s / TR).
        bins (int): number of histogram bins for the FCD distribution.
        get_FCD_matrix (bool): also return the full FCD matrix.
    Returns:
        ``fcd_mat`` (Tensor): ``[M, window_num, window_num]``.
        ``fcd_hist`` (Tensor): ``[bins, M]`` histogram counts.
    """
    set_torch_default()

    if bold.ndim == 2:
        bold = bold.unsqueeze(1)

    N = bold.shape[0]
    M = bold.shape[1]
    t_len = bold.shape[2]
    window_num = t_len - window_size + 1
    if t_len < window_size:
        raise Exception("The length of bold signal is shorter than the window size!")

    # Sliding-window FCs.
    fc_list = torch.zeros(M, window_num, N, N)
    for t in range(0, window_num):
        bold_single = bold[:, :, t:t + window_size]
        fc_list[:, t, :, :] = calc_FC(bold_single)

    # FCD = correlation between sliding-window FCs (vectorised upper triangle).
    fcd_mat = torch.zeros(M, window_num, window_num)
    for m in range(M):
        fc_vectors = vectorize_fc(fc_list[m])  # [window_num, N*(N-1)/2]
        fcd_mat[m] = torch.corrcoef(fc_vectors)

    # Off-diagonal entries -> histogram.
    fcd_mask = torch.ones(window_num, window_num, dtype=torch.bool)
    fcd_mask = torch.triu(fcd_mask, 1)
    fcd_vec = fcd_mat[:, fcd_mask]
    fcd_hist = torch.ones(bins, M)
    for hist_i in range(M):
        fcd_hist[:, hist_i] = torch.histc(
            fcd_vec[hist_i], bins=bins, min=-1.0, max=1.0
        )

    if fcd_hist.shape[1] == 1:
        fcd_hist = fcd_hist.squeeze(1)
    if fcd_mat.shape[0] == 1:
        fcd_mat = fcd_mat.squeeze(0)

    # Note: `get_FCD_matrix` is honoured by callers via the second return value.
    return fcd_mat, fcd_hist


def fcd_hist_to_cdf(fcd_hist, normalize=True):
    """Cumulative-sum convert an FCD histogram into an FCD-CDF."""
    fcd_cdf = torch.cumsum(fcd_hist, dim=0)
    if normalize:
        if fcd_hist.ndim == 1:
            fcd_cdf = fcd_cdf / fcd_cdf[-1]
        else:
            fcd_cdf = fcd_cdf / fcd_cdf[-1:, :]
    return fcd_cdf


def get_fcd_max_count(TR, scan_length):
    """Number of distinct sliding-window pairs at the standard 60-s window."""
    num_time_points = round(scan_length / TR)
    window_size = round(60 / TR)
    num_windows = num_time_points - window_size + 1
    max_count = num_windows * (num_windows - 1) // 2
    return max_count


# --------------------------------------------------------------------------- #
# Cost components
# --------------------------------------------------------------------------- #
def calc_FC_losses(fc_sim, emp_fc):
    """FC correlation (``1 - r``) and FC L1 costs.

    Args:
        fc_sim (Tensor): ``[M, N, N]`` or ``[N, N]`` simulated FC matrices.
        emp_fc (Tensor): ``[N, N]`` empirical FC matrix.
    Returns:
        dict with keys ``corr_loss``, ``l1_loss`` (= ``|mean(FC_sim)-mean(FC_emp)|``),
        ``old_l1_loss``, ``MAE_l1_loss``.
    """
    if fc_sim.ndim == 2:
        fc_sim = fc_sim.unsqueeze(0)
    M = fc_sim.shape[0]

    vec_emp = vectorize_fc(emp_fc)
    vec_emp = vec_emp.unsqueeze(0).expand(M, -1)
    vec_sim = vectorize_fc(fc_sim)

    # "L1" in the paper = absolute difference of FC means (NOT MAE).
    old_l1_loss = torch.abs(
        torch.mean(vec_emp, dim=1) - torch.mean(vec_sim, dim=1)
    )
    MAE_l1_loss = torch.mean(torch.abs(vec_emp - vec_sim), dim=1)

    vec_3d = torch.zeros(M, 2, vec_emp.shape[1])
    vec_3d[:, 0, :] = vec_sim
    vec_3d[:, 1, :] = vec_emp
    corr = torch_corr_3D(vec_3d)
    corr_loss = 1 - corr

    return {
        "corr_loss": corr_loss,
        "l1_loss": old_l1_loss,
        "old_l1_loss": old_l1_loss,
        "MAE_l1_loss": MAE_l1_loss,
    }


def calc_FCD_losses(sim_fcd_hist, emp_fcd_cdf):
    """KS distance between simulated FCD histograms and an empirical FCD CDF.

    Args:
        sim_fcd_hist (Tensor): ``[bins, M]`` simulated FCD histograms.
        emp_fcd_cdf  (Tensor): ``[bins, 1]`` empirical FCD CDF (already normalised).
    Returns:
        dict ``{"ks_loss": Tensor[M]}``.
    """
    M = sim_fcd_hist.shape[1]
    sim_fcd_cum = fcd_hist_to_cdf(sim_fcd_hist, normalize=True)
    emp_fcd_cdf_expand = emp_fcd_cdf.expand(-1, M)
    ks_dif = torch.abs(sim_fcd_cum - emp_fcd_cdf_expand)
    ks_loss = torch.max(ks_dif, dim=0)[0]
    return {"ks_loss": ks_loss}


def calc_all_loss_from_fc_fcd(fc_sim, fcd_hist_sim, emp_fc, emp_fcd_cdf):
    """Combine FC and FCD losses into a single dict (one value per param set)."""
    losses = calc_FC_losses(fc_sim, emp_fc)
    losses.update(calc_FCD_losses(fcd_hist_sim, emp_fcd_cdf))
    return losses
