"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch
import scipy.io as spio
import numpy as np
import pandas as pd
import os


def read_predictor_constants(config):
    """Read the cost predictor settings that DELSSOME mode needs from a CMA-ES config.

    Args:
        config (configparser.ConfigParser): the parsed .ini config

    Returns:
        tuple: (predictor_save_path (str), d_transformer (int))

    Raises:
        Exception: if the config has no [Deep learning Model Constants] section
    """
    if 'Deep learning Model Constants' not in config:
        raise Exception(
            'DELSSOME mode needs a [Deep learning Model Constants] section in the config, holding '
            '"predictor_save_path" and "d_transformer". See '
            'examples/group_level/config/example_CMAES_pFIC_DELSSOME.ini for a working example, or '
            'pass --train_mode Euler to score candidates by direct simulation instead.')
    dl_model_constants = config['Deep learning Model Constants']
    return dl_model_constants['predictor_save_path'], int(dl_model_constants['d_transformer'])


def convert_csv_to_tensor(input_dir, file_name):
    """Convert the contents of csv file to torch.Tensor.
    Check the existence of file and raise Exception if file not exist.

    Args:
        input_dir (str): input directory
        file_name (str): csv file name

    Returns:
        Tensor: converted tensor
    """
    if not file_name.endswith('.csv'):
        raise Exception("File name should contain .csv")
    file_path = os.path.join(input_dir, file_name)
    if not os.path.exists(file_path):
        raise Exception(f"{file_name} not in the input directory.")
    file_content = np.array(pd.read_csv(file_path, header=None, index_col=False))
    file_content = torch.as_tensor(file_content)
    return file_content


def convert_mat_to_tensor(input_dir, file_name, mat_key):
    """Convert the contents of mat file to torch.Tensor.
    Check the existence of file and raise Exception if file not exist.

    Args:
        input_dir (str): input directory
        file_name (str): mat file name

    Returns:
        Tensor: converted tensor
    """
    if not file_name.endswith('.mat'):
        raise Exception("File name should contain .mat")
    file_path = os.path.join(input_dir, file_name)
    if not os.path.exists(file_path):
        raise Exception(f"{file_name} not in the input directory.")
    file_content = spio.loadmat(file_path)
    if mat_key not in file_content:
        raise Exception("Key does not exist.")
    file_content = torch.as_tensor(file_content[mat_key])
    return file_content


def set_torch_default(GPU_index):
    """Set default float dtype to float64 and print the active CUDA device (if available)."""
    if torch.cuda.is_available():
        print('CUDA Device: ', torch.cuda.get_device_name())
    torch.set_default_dtype(torch.float64)


def CBIG_corr(s_series, t_series=None):
    """Compute the Pearson correlation for s_series (self_correlation) or between s and t.
    [IMPORTANT] features are in dim=0 (intrinsically vertical in matrix) and labels are in dim=1

    Args:
        s_series (tensor): [features, labels]
        t_series (tensor, optional): [features, t_labels]. Defaults to None.

    Returns:
        array: [labels, labels] or [labels, t_labels]
    """
    s_series = s_series - torch.mean(s_series, dim=0)
    s_series = s_series / torch.linalg.norm(s_series, dim=0)
    if t_series is None:
        return torch.matmul(torch.transpose(s_series, 0, 1), s_series)
    else:
        t_series = t_series - torch.mean(t_series, dim=0)
        t_series = t_series / torch.linalg.norm(t_series, dim=0)
        return torch.matmul(torch.transpose(s_series, 0, 1), t_series)


def parameterize_myelin_rsfc(myelin, rsfc_gradient, param_10):
    """get full parameters from parameterization coefficients

    Args:
        myelin (tensor): myelin profile of [N, 1]
        rsfc_gradient (tensor): RSFC gradient profile of [N, 1]
        param_10 (tensor): parametrization coefficients, [10, param_sets]

    Returns:
        tensor: parameters for MFM model, [3N+1, param_sets]
    """
    w_EE = param_10[0] + param_10[1] * myelin + param_10[2] * rsfc_gradient
    w_EI = param_10[3] + param_10[4] * myelin + param_10[5] * rsfc_gradient
    G = param_10[6]
    sigma = param_10[7] + param_10[8] * myelin + param_10[9] * rsfc_gradient
    return torch.vstack((w_EE, w_EI, G, sigma))


def parameterize_myelin_rsfc_MFM2013(myelin, rsfc_gradient, param_10):
    """get full parameters from parameterization coefficients

    Args:
        myelin (tensor): myelin profile of [N, 1]
        rsfc_gradient (tensor): RSFC gradient profile of [N, 1]
        param_10 (tensor): parametrization coefficients, [10, param_sets]

    Returns:
        tensor: parameters for MFM model, [3N+1, param_sets]
    """
    w = param_10[0] + param_10[1] * myelin + param_10[2] * rsfc_gradient
    I_0 = param_10[3] + param_10[4] * myelin + param_10[5] * rsfc_gradient
    G = param_10[6]
    sigma = param_10[7] + param_10[8] * myelin + param_10[9] * rsfc_gradient
    return torch.vstack((w, I_0, G, sigma))


def parameterize_myelin_rsfc_Hopf(myelin, rsfc_gradient, param_10):
    """get full parameters from parameterization coefficients

    Args:
        myelin (tensor): myelin profile of [N, 1]
        rsfc_gradient (tensor): RSFC gradient profile of [N, 1]
        param_10 (tensor): parametrization coefficients, [10, param_sets]

    Returns:
        tensor: parameters for MFM model, [3N+1, param_sets]
    """
    a = param_10[0] + param_10[1] * myelin + param_10[2] * rsfc_gradient
    omega = param_10[3] + param_10[4] * myelin + param_10[5] * rsfc_gradient
    G = param_10[6]
    sigma = param_10[7] + param_10[8] * myelin + param_10[9] * rsfc_gradient
    return torch.vstack((a, omega, G, sigma))


def torch_corr_3D(vec_3d):
    """
    Compute the correlation coefficient parallel for 3D vector.
    And the 2nd dimension must be 2, standing for X and Y.
    :param vec_3d: [M, 2, len]
    :return: The corr_coef for X and Y. corr_coef: [M,]
    """
    ex_ey = torch.mean(vec_3d, dim=-1)
    std = torch.std(vec_3d, dim=-1, unbiased=False)
    xy = vec_3d[:, 0, :] * vec_3d[:, 1, :]
    e_xy = torch.mean(xy, dim=-1)
    cov = e_xy - ex_ey[:, 0] * ex_ey[:, 1]
    corr = cov / (std[:, 0] * std[:, 1])
    return corr


def all_loss_calculate_from_fc_fcd(fc_sim, fcd_hist_sim, fc_emp, fcd_cdf_emp):
    """
    Calculate corr_loss, L1_loss and KS_loss from simulated FC matrix and FCD histogram
    :param fc_sim: [M, N, N], sets of N x N simulated FC matrices
    :param fcd_hist: [fcd_hist_bins, M], no need to cumsum or normalize. Will do automatically in KS_cost
    :param fc_emp: [N, N], the empirical FC matrix
    :param fcd_cdf_emp: [fcd_hist_bins, 1], the empirical FCD CDF
    :return: [M], sets of loss
    """
    corr, L1_loss = FC_correlation_n_L1_cost(fc_sim, fc_emp)
    corr_loss = 1 - corr
    ks_loss = KS_cost(fcd_hist_sim, fcd_cdf_emp)
    total_loss = corr_loss + L1_loss + ks_loss
    return total_loss, corr_loss, L1_loss, ks_loss


def all_loss_calculate_from_bold(bold, fc_emp, fcd_cdf_emp):
    """
    Calculate corr_loss, L1_loss and KS_loss from BOLD signals
    :param bold: [N, M, len], the simulated BOLD signals
    :param fc_emp: [N, N], the empirical FC matrix
    :param fcd_cdf_emp: [fcd_cdf_emp, 1]. Has been done cumulative summation and normalization
    (dividing by emp_fcd_cum[-1:, :])
    :return: total_loss: [M]
    """
    fc_sim = FC_calculate(bold)
    corr, L1_loss = FC_correlation_n_L1_cost(fc_sim, fc_emp)
    corr_loss = 1 - corr
    _, fcd_hist = FCD_calculate(bold)
    ks_loss = KS_cost(fcd_hist, fcd_cdf_emp)
    total_loss = corr_loss + L1_loss + ks_loss
    return total_loss


def FC_calculate(bold):
    """
    Calculate FC matrix for M sets of BOLD signals
    :param bold: [N, M, t_len]
    :return: FC matrix [M, N, N]
    """
    N = bold.shape[0]
    M = bold.shape[1]
    fc_mat = torch.zeros((M, N, N))
    for i in range(M):
        fc_mat[i] = torch.corrcoef(bold[:, i, :])
    return fc_mat


def FC_correlation_n_L1_cost(fc_sim, fc_emp):
    """
    Compute the FC correlation and L1 cost for all M sets.
    [ATTENTION] here L1 is not the real L1 (absolute difference of values and get the mean),
    but the absolute difference of two mean values.
    :param fc_sim: [M, N, N]
    :param fc_emp: [N, N]
    :return: corr: [M,]; L1: [M,]
    """
    M = fc_sim.shape[0]
    N = fc_sim.shape[1]

    # Extract the upper triangular part
    mask = torch.ones(N, N, dtype=torch.bool)
    mask = torch.triu(mask, 1)
    vec_emp = fc_emp[mask]
    vec_emp = vec_emp.unsqueeze(0).expand(M, -1)  # [M, len]
    vec_sim = fc_sim[:, mask]

    # L1 version 1: abs(mean)
    L1_cost = torch.abs(torch.mean(vec_emp, dim=1) - torch.mean(vec_sim, dim=1))

    vec_3d = torch.zeros(M, 2, vec_emp.shape[1])
    vec_3d[:, 0, :] = vec_sim
    vec_3d[:, 1, :] = vec_emp

    corr = torch_corr_3D(vec_3d)
    return corr, L1_cost


def FC_correlation_single(fc_1, fc_2):
    """
    Compute the correlation between two FC matrix.
    By flatten their upright part to two vectors and use Pearson's correlation
    :param fc_1: [N, N]
    :param fc_2: [N, N]
    :return:
    """
    N = fc_1.shape[0]
    mask = np.ones((N, N)).astype(bool)
    mask = np.triu(mask, 1)
    vec = np.zeros((2, (N * N - N) // 2))
    vec[0] = fc_1[mask]
    vec[1] = fc_2[mask]
    cor_fc = np.corrcoef(vec)
    return cor_fc[0, 1]


def FCD_calculate(bold, window_size=83, bins=10000):
    """
    Moving windows and calculating the FC matrix for every window_size,
    then calculating the correlation between these FC matrix.
    :param bold: [N, M, t_len]
    :param window_size: The length of sliding window
    :param bins: The histogram bins
    :return: FCD matrix [M, window_num, window_num]. window_num = t_len - window_size + 1
            FCD histogram: [bins, M]
    """

    N = bold.shape[0]
    M = bold.shape[1]
    t_len = bold.shape[2]
    window_num = t_len - window_size + 1
    if t_len < window_size:
        raise Exception("The length of bold signal is shorter than the window size!")

    # Compute sliding window FC
    fc_list = torch.zeros(M, window_num, N, N)
    for t in range(0, window_num):
        bold_single = bold[:, :, t:t + window_size]
        fc_list[:, t, :, :] = FC_calculate(bold_single)

    # Compute pairwise correlation between FCs
    fcd_mat = torch.zeros(M, window_num, window_num)
    fc_mask = torch.ones(N, N, dtype=torch.bool)
    fc_mask = torch.triu(fc_mask, 1)
    for m in range(M):
        fcd_mat[m] = torch.corrcoef(fc_list[m, :, fc_mask])

    # Calculate the FCD histogram
    fcd_mask = torch.ones(window_num, window_num, dtype=torch.bool)
    fcd_mask = torch.triu(fcd_mask, 1)
    fcd_vec = fcd_mat[:, fcd_mask]
    fcd_hist = torch.ones(bins, M)
    for hist_i in range(M):
        fcd_hist[:, hist_i] = torch.histc(fcd_vec[hist_i], bins=bins, min=-1., max=1.)

    return fcd_mat, fcd_hist


def KS_cost(fcd_hist_sim, fcd_cdf_emp):
    """
    Calculate the KS cost between the simulated FCD histogram and the empirical one
    :param sim_fcd_hist: [bins, M]
    :param emp_fcd_cum: [bins, 1]. Has been done cumulative summation and normalization
    (divided by emp_fcd_cum[-1:, :])
    :return: KS_cost: [M]
    """
    M = fcd_hist_sim.shape[1]
    sim_fcd_cum = torch.cumsum(fcd_hist_sim, dim=0)
    sim_fcd_cum = sim_fcd_cum / sim_fcd_cum[-1:, :]
    emp_fcd_cum_expand = fcd_cdf_emp.expand(-1, M)
    ks_dif = torch.abs(sim_fcd_cum - emp_fcd_cum_expand)
    ks_cost = torch.max(ks_dif, dim=0)[0]
    return ks_cost
