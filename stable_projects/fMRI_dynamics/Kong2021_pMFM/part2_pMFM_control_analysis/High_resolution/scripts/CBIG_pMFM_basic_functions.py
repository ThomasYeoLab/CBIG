# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import csv
import math
import time
import numpy as np
import torch
import scipy.io as sio
'''**********  Functions for computing simulated BOLD signals   ************'''


def CBIG_mfm_multi_simulation(parameter, sc_mat, t_epochlong, n_dup):
    '''
    Function used to generate the simulated BOLD signal using mean field model
    and hemodynamic model
    Each parameter set is ussed to simulated multiple times to get stable
    result
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    Each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        sc_mat:     N*N structural connectivity matrix
        t_epochlong:total simulated time
        n_dup:      Number of times each parameter set is simulated
    Returns:
        bold_d:     simulated BOLD signal
    '''

    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    # Initializing system parameters
    kstart = 0.
    t_pre = 60 * 2
    kend = t_pre + 60 * t_epochlong
    d_t = 0.01
    t_bold = 0.72

    # Setting sampling ratio
    k_p = torch.arange(kstart, kend + d_t, d_t)
    n_num = parameter.shape[1]
    n_set = n_dup * n_num
    parameter = parameter.repeat(1, n_dup)
    n_nodes = sc_mat.shape[0]
    n_samples = k_p.shape[0]

    # Initializing neural activity
    y_t = torch.zeros((n_nodes, n_set))
    d_y = torch.zeros((n_nodes, n_set))

    # Initializing hemodynamic activity
    f_mat = torch.ones((n_nodes, n_set, 4))
    z_t = torch.zeros((n_nodes, n_set))
    f_t = torch.ones((n_nodes, n_set))
    v_t = torch.ones((n_nodes, n_set))
    q_t = torch.ones((n_nodes, n_set))
    f_mat[:, :, 0] = z_t
    y_t[:, :] = 0.001

    # Wiener process
    w_coef = parameter[2 * n_nodes + 1:3 * n_nodes + 1, :] / math.sqrt(0.001)
    if w_coef.shape[0] == 1:
        w_coef = w_coef.repeat(n_nodes, 1)
    w_l = k_p.shape[0]
    d_w = math.sqrt(d_t) * torch.randn(n_dup, n_nodes, w_l + 1000)
    p_costant = 0.34
    v_0 = 0.02
    k_1 = 4.3 * 28.265 * 3 * 0.0331 * p_costant
    k_2 = 0.47 * 110 * 0.0331 * p_costant
    k_3 = 0.53
    count = 0
    y_bold = torch.zeros((n_nodes, n_set, int(n_samples / (t_bold / d_t) + 1)))

    # Warm up
    start = time.time()
    for i in range(1000):
        d_y = CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        noise_level = d_w[:, :, i].repeat(1, 1, n_num).contiguous().view(
            -1, n_nodes)
        y_t = y_t + d_y * d_t + w_coef * torch.transpose(noise_level, 0, 1)

    # Main body: calculation
    for i in range(n_samples):
        d_y = CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        noise_level = d_w[:, :, i + 1000].\
            repeat(1, 1, n_num).contiguous().view(-1, n_nodes)
        y_t = y_t + d_y * d_t + w_coef * torch.transpose(noise_level, 0, 1)
        d_f = CBIG_mfm_rfMRI_BW_ode(y_t, f_mat)
        f_mat = f_mat + d_f * d_t
        z_t, f_t, v_t, q_t = torch.chunk(f_mat, 4, dim=2)
        y_bold_temp = 100 / p_costant * v_0 * \
            (k_1 * (1 - q_t) + k_2 * (1 - q_t / v_t) + k_3 * (1 - v_t))
        y_bold[:, :, count] = y_bold_temp[:, :, 0]
        count = count + ((i + 1) % (t_bold / d_t) == 0) * 1
    elapsed = time.time() - start
    print('The time used for calculating simulated BOLD signal is: ', elapsed)

    # Downsampling
    cut_index = int(t_pre / t_bold)
    bold_d = y_bold[:, :, cut_index + 1:y_bold.shape[2]]

    return bold_d


def CBIG_mfm_single_simulation(parameter, sc_mat, t_epochlong, d_t=0.01):
    '''
    Function used to generate the simulated BOLD signal using mean field model
    and hemodynamic model
    Each parameter set is ussed to simulated one time
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets
                    Each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        sc_mat:     N*N structural connectivity matrix
        t_epochlong:total simulated time
    Returns:
        bold_d:     simulated BOLD signal
    '''

    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    # Initializing system parameters
    kstart = 0.
    t_pre = 60 * 2
    kend = t_pre + 60 * t_epochlong
    t_bold = 0.72

    # sampling ratio
    k_p = torch.arange(kstart, kend + d_t, d_t)
    n_nodes = sc_mat.shape[0]
    n_samples = k_p.shape[0]
    n_set = parameter.shape[1]

    # Initializing neural activity
    y_t = torch.zeros((n_nodes, n_set))
    d_y = torch.zeros((n_nodes, n_set))

    # Initializing hemodynamic activity
    f_mat = torch.ones((n_nodes, n_set, 4))
    z_t = torch.zeros((n_nodes, n_set))
    f_t = torch.ones((n_nodes, n_set))
    v_t = torch.ones((n_nodes, n_set))
    q_t = torch.ones((n_nodes, n_set))
    f_mat[:, :, 0] = z_t
    y_t[:, :] = 0.001

    # Wiener process
    w_coef = parameter[2 * n_nodes + 1:3 * n_nodes + 1, :] / math.sqrt(0.001)
    if w_coef.shape[0] == 1:
        w_coef = w_coef.repeat(n_nodes, 1)
    p_costant = 0.34
    v_0 = 0.02
    k_1 = 4.3 * 28.265 * 3 * 0.0331 * p_costant
    k_2 = 0.47 * 110 * 0.0331 * p_costant
    k_3 = 0.53
    count = 0
    y_bold = torch.zeros((n_nodes, n_set, int(n_samples / (t_bold / d_t) + 1)))

    # Warm up
    start = time.time()
    for i in range(1000):
        d_y = CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        y_t = y_t + d_y * d_t + w_coef * \
            torch.randn(n_nodes, n_set) * math.sqrt(d_t)

    # Main body: calculation
    for i in range(n_samples):
        d_y = CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        random_num = torch.randn(n_nodes, n_set)
        y_t = y_t + d_y * d_t + w_coef * random_num * math.sqrt(d_t)
        d_f = CBIG_mfm_rfMRI_BW_ode(y_t, f_mat)
        f_mat = f_mat + d_f * d_t
        z_t, f_t, v_t, q_t = torch.chunk(f_mat, 4, dim=2)
        y_bold_temp = 100 / p_costant * v_0 * \
            (k_1 * (1 - q_t) + k_2 * (1 - q_t / v_t) + k_3 * (1 - v_t))
        y_bold[:, :, count] = y_bold_temp[:, :, 0]
        count = count + ((i + 1) % (t_bold / d_t) == 0) * 1
    elapsed = time.time() - start
    print('The time used for calculating simulated BOLD signal is: ', elapsed)

    # Downsampling
    cut_index = int(t_pre / t_bold)
    bold_d = y_bold[:, :, cut_index + 1:y_bold.shape[2]]

    return bold_d


def CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat):
    '''
    This function is to calcualte the derivatives of synaptic gating variable S
    Args:
        y_t:        N*M matrix represents synaptic gating variable
                    N is the number of ROI
                    M is the number of candidate parameter sets
        parameter:  (N*3+1)*M matrix.
                    Each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        sc_mat:     N*N structural connectivity matrix
    Returns:
        dy:         N*M matrix represents derivatives of synaptic gating
                    variable S
    '''

    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    # Parameters for inputs and couplings
    number_roi = sc_mat.shape[0]
    J = 0.2609
    w = parameter[0:number_roi, :]
    G = parameter[2 * number_roi, :]
    I0 = parameter[number_roi:2 * number_roi, :]

    # Parameters for firing rate
    a = 270
    b = 108
    d = 0.154

    # Parameters for synaptic activity/currents
    tau_s = 0.1
    gamma_s = 0.641

    # Total input x
    x = J * w * y_t + J * G.repeat(number_roi, 1) * torch.mm(sc_mat, y_t) + I0

    # Population firing rate
    H = (a * x - b) / (1 - torch.exp(-d * (a * x - b)))

    # Synaptic activity
    dy = -1 / tau_s * y_t + gamma_s * (1 - y_t) * H

    return dy


def CBIG_mfm_rfMRI_BW_ode(y_t, F):
    '''
    This fucntion is to implement the hemodynamic model
    Args:
        y_t:        N*M matrix represents synaptic gating variable
                    N is the number of ROI
                    M is the number of candidate parameter sets
        F:          Hemodynamic activity variables
    Returns:
        dF:         Derivatives of hemodynamic activity variables
    '''

    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    # Hemodynamic model parameters
    beta = 0.65
    gamma = 0.41
    tau = 0.98
    alpha = 0.33
    p_constant = 0.34
    n_nodes = y_t.shape[0]
    n_set = y_t.shape[1]

    # Calculate derivatives
    dF = torch.zeros((n_nodes, n_set, 4))
    dF[:, :, 0] = y_t - beta * F[:, :, 0] - gamma * (F[:, :, 1] - 1)
    dF[:, :, 1] = F[:, :, 0]
    dF[:, :, 2] = 1 / tau * (F[:, :, 1] - F[:, :, 2]**(1 / alpha))
    dF[:, :, 3] = 1 / tau * (F[:, :, 1] / p_constant * (1 - (1 - p_constant)**
                                                        (1 / F[:, :, 1])) -
                             F[:, :, 3] / F[:, :, 2] * F[:, :, 2]**(1 / alpha))
    return dF


def CBIG_FCcorrelation_multi_simulation(emp_fc, bold_d, n_dup):
    '''
    This function is to calculate the FC correlation cost for multiple
    simulation BOLD signal results
    Args:
        emp_fc:     N*N group level FC matrix
                    N is number of ROI
        bold_d:     simulated BOLD signal
        n_dup:      Number of times each parameter set is simulated
    Returns:
        corr_cost:   FC correlation cost
    '''

    fc_timestart = time.time()

    # Calculate vectored simulated FC
    n_set = bold_d.shape[1]
    n_num = int(n_set / n_dup)
    n_nodes = emp_fc.shape[0]
    fc_mask = torch.triu(torch.ones(n_nodes, n_nodes), 1) == 1
    vect_len = int(n_nodes * (n_nodes - 1) / 2)
    sim_fc_vector = torch.zeros(n_set, vect_len)
    for i in range(n_set):
        sim_fc = torch_corr(bold_d[:, i, :])
        sim_fc_vector[i, :] = sim_fc[fc_mask]

    # Average the simulated FCs with same parameter set
    sim_fc_vector[sim_fc_vector != sim_fc_vector] = 0
    sim_fc_num = torch.zeros(n_num, vect_len)
    sim_fc_den = torch.zeros(n_num, 1)
    for k in range(n_dup):
        sim_fc_num = sim_fc_num + sim_fc_vector[k * n_num:(k + 1) * n_num, :]
        sim_fc_den = sim_fc_den + \
            (sim_fc_vector[k * n_num:(k + 1) * n_num, 0:1] != 0).float()
    sim_fc_den[sim_fc_den == 0] = np.nan
    sim_fc_ave = sim_fc_num / sim_fc_den

    # Calculate FC correlation
    emp_fcm = emp_fc[fc_mask].repeat(n_num, 1)
    corr_mass = torch_corr2(torch_arctanh(sim_fc_ave), torch_arctanh(emp_fcm))
    corr_cost = torch.diag(corr_mass)
    corr_cost = corr_cost.cpu().numpy()
    corr_cost = 1 - corr_cost
    corr_cost[np.isnan(corr_cost)] = 10

    fc_elapsed = time.time() - fc_timestart
    print('Time using for calcualting FC correlation cost: ', fc_elapsed)

    return corr_cost


def CBIG_FCcorrelation_single_simulation(emp_fc, bold_d, n_dup):
    '''
    This function is to calculate the FC correlation cost for single simulation
    BOLD signal result
    Args:
        emp_fc:     N*N group level FC matrix
                    N is number of ROI
        bold_d:     simulated BOLD signal
    Returns:
        corr_cost:   FC correlation cost
    '''

    fc_timestart = time.time()

    # Calculate vectored simulated FC
    n_set = bold_d.shape[1]
    n_nodes = emp_fc.shape[0]
    fc_mask = torch.triu(torch.ones(n_nodes, n_nodes), 1) == 1
    vect_len = int(n_nodes * (n_nodes - 1) / 2)
    sim_fc_vector = torch.zeros(n_set, vect_len)
    for i in range(n_set):
        sim_fc = torch_corr(bold_d[:, i, :])
        sim_fc_vector[i, :] = sim_fc[fc_mask]

    # Calculate FC correlation
    sim_fc_numpy = sim_fc_vector.cpu().numpy()
    emp_fc_numpy = emp_fc[fc_mask].cpu().numpy()
    time_dup = int(n_set / n_dup)
    corr_cost = np.zeros(time_dup)

    for t in range(time_dup):
        sim_fc_numpy_temp = sim_fc_numpy[t * n_dup:(t + 1) * n_dup, :]
        sim_fc_mean = np.nanmean(sim_fc_numpy_temp, 0)
        corrmean_temp = np.corrcoef(
            np.arctanh(sim_fc_mean), np.arctanh(emp_fc_numpy))
        corr_cost[t] = 1 - corrmean_temp[1, 0]

    fc_elapsed = time.time() - fc_timestart
    print('Time using for calcualting FC correlation cost: ', fc_elapsed)
    return corr_cost


def CBIG_FCDKSstat_multi_simulation(emp_ks, bold_d, n_dup):
    '''
    This function is to calculate the FCD KS statistics cost for multiple
    simulation BOLD signal results
    Args:
        emp_ks:     Group level KS statistics for empirical data
        bold_d:     simulated BOLD signal
        n_dup:      Number of times each parameter set is simulated
    Returns:
        ks_cost:   FCD KS statistics cost
    '''

    fcd_timestart = time.time()

    # Initializing the FC and FCD masks
    n_set = bold_d.shape[1]
    n_num = int(n_set / n_dup)
    n_nodes = bold_d.shape[0]
    window_size = 83
    time_lengh = 1200 - window_size + 1
    sub_num = 10
    resid_num = n_set % sub_num
    fc_edgenum = int(n_nodes * (n_nodes - 1) / 2)
    fc_mask = torch.triu(torch.ones(n_nodes, n_nodes), 1) == 1
    fc_maskm = torch.zeros(n_nodes * sub_num,
                           n_nodes * sub_num).type(torch.cuda.ByteTensor)
    for i in range(sub_num):
        fc_maskm[n_nodes * i:n_nodes * (i + 1), n_nodes * i:n_nodes *
                 (i + 1)] = fc_mask

    fc_mask_resid = torch.zeros(n_nodes * resid_num, n_nodes * resid_num).type(
        torch.cuda.ByteTensor)
    for i in range(resid_num):
        fc_mask_resid[n_nodes * i:n_nodes * (i + 1), n_nodes * i:n_nodes *
                      (i + 1)] = fc_mask

    fcd_mask = torch.triu(torch.ones(time_lengh, time_lengh), 1) == 1

    # Calculating CDF for simualted FCD matrices
    fcd_hist = torch.ones(10000, n_set).cpu()
    fc_mat = torch.zeros(fc_edgenum, sub_num, time_lengh)
    batch_num = math.floor(n_set / sub_num)
    fc_resid = torch.zeros(fc_edgenum, resid_num, time_lengh)

    for b in range(batch_num):
        bold_temp = bold_d[:, b * sub_num:(b + 1) * sub_num, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, 1200)
        for i in range(0, time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_maskm]
            fc_mat[:, :, i] = torch.transpose(
                cor_temp.view(sub_num, fc_edgenum), 0, 1)

        for j in range(0, sub_num):
            fcd_temp = torch_corr(torch.transpose(fc_mat[:, j, :], 0, 1))
            fcd_hist[:, j + b * sub_num] = torch.histc(
                fcd_temp[fcd_mask].cpu(), 10000, (0.0001 - 1), 1)

    if resid_num != 0:
        bold_temp = bold_d[:, batch_num * sub_num:n_set, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, 1200)
        for i in range(time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_mask_resid]
            fc_resid[:, :, i] = torch.transpose(
                cor_temp.view(resid_num, fc_edgenum), 0, 1)

        for j in range(resid_num):
            fcd_temp = torch_corr(torch.transpose(fc_resid[:, j, :], 0, 1))
            fcd_hist[:, j + sub_num * batch_num] = torch.histc(
                fcd_temp[fcd_mask].cpu(), 10000, (0.0001 - 1), 1)

    fcd_histcum = np.cumsum(fcd_hist.numpy(), 0)
    fcd_histcumM = fcd_histcum.copy()
    fcd_histcumM[:, fcd_histcum[-1, :] != emp_ks[-1, 0]] = 0

    # Calculating KS statistics cost
    fcd_histcum_temp = np.zeros((10000, n_num))
    fcd_histcum_num = np.zeros((1, n_num))
    for k in range(n_dup):
        fcd_histcum_temp = fcd_histcum_temp + \
            fcd_histcumM[:, k * n_num:(k + 1) * n_num]
        fcd_histcum_num = fcd_histcum_num + \
            (fcd_histcumM[-1, k * n_num:(k + 1) * n_num] == emp_ks[-1, 0])
    fcd_histcum_ave = fcd_histcum_temp / fcd_histcum_num
    ks_diff = np.abs(fcd_histcum_ave - np.tile(emp_ks, [1, n_num]))
    ks_cost = ks_diff.max(0) / emp_ks[-1, 0]
    ks_cost[fcd_histcum_ave[-1, :] != emp_ks[-1, 0]] = 10

    fcd_elapsed = time.time() - fcd_timestart
    print('Time using for calcualting FCD KS statistics cost: ', fcd_elapsed)
    return ks_cost


def CBIG_FCDKSstat_single_simulation(emp_ks, bold_d, n_dup, window_size=83):
    '''
    This function is to calculate the FCD KS statistics cost for single
    simulation BOLD signal results
    Args:
        emp_ks:     Group level KS statistics for empirical data
        bold_d:     simulated BOLD signal
    Returns:
        ks_cost:    FCD KS statistics cost
    '''

    fcd_timestart = time.time()

    # Initializing the FC and FCD masks
    n_set = bold_d.shape[1]
    n_nodes = bold_d.shape[0]
    time_lengh = 1200 - window_size + 1
    sub_num = 10
    resid_num = n_set % sub_num
    fc_edgenum = int(n_nodes * (n_nodes - 1) / 2)
    fc_mask = torch.triu(torch.ones(n_nodes, n_nodes), 1) == 1
    fc_maskm = torch.zeros(n_nodes * sub_num,
                           n_nodes * sub_num).type(torch.cuda.ByteTensor)

    for i in range(sub_num):
        fc_maskm[n_nodes * i:n_nodes * (i + 1), n_nodes * i:n_nodes *
                 (i + 1)] = fc_mask

    fc_mask_resid = torch.zeros(n_nodes * resid_num, n_nodes * resid_num).type(
        torch.cuda.ByteTensor)
    for i in range(resid_num):
        fc_mask_resid[n_nodes * i:n_nodes * (i + 1), n_nodes * i:n_nodes *
                      (i + 1)] = fc_mask

    fcd_mask = torch.triu(torch.ones(time_lengh, time_lengh), 1) == 1

    # Calculating CDF for simualted FCD matrices
    fcd_hist = torch.ones(10000, n_set).cpu()
    fc_mat = torch.zeros(fc_edgenum, sub_num, time_lengh)
    batch_num = int(n_set / sub_num)
    fc_resid = torch.zeros(fc_edgenum, resid_num, time_lengh)

    for b in range(batch_num):
        bold_temp = bold_d[:, b * sub_num:(b + 1) * sub_num, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, 1200)
        for i in range(0, time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_maskm]
            fc_mat[:, :, i] = torch.transpose(
                cor_temp.view(sub_num, fc_edgenum), 0, 1)

        for j in range(0, sub_num):
            fcd_temp = torch_corr(torch.transpose(fc_mat[:, j, :], 0, 1))
            fcd_hist[:, j + b * sub_num] = torch.histc(
                fcd_temp[fcd_mask].cpu(), 10000, (0.0001 - 1), 1)

    if resid_num != 0:
        bold_temp = bold_d[:, batch_num * sub_num:n_set, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, 1200)
        for i in range(time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_mask_resid]
            fc_resid[:, :, i] = torch.transpose(
                cor_temp.view(resid_num, fc_edgenum), 0, 1)

        for j in range(resid_num):
            fcd_temp = torch_corr(torch.transpose(fc_resid[:, j, :], 0, 1))
            fcd_hist[:, j + sub_num * batch_num] = torch.histc(
                fcd_temp[fcd_mask].cpu(), 10000, (0.0001 - 1), 1)

    fcd_histcum = np.cumsum(fcd_hist.numpy(), 0)

    # Calculating KS statistics cost
    time_dup = int(n_set / n_dup)
    ks_cost = np.zeros(time_dup)
    for t in range(time_dup):
        fcd_hist_temp = fcd_histcum[:, t * n_dup:(t + 1) * n_dup]
        fcd_histcum_nn = fcd_hist_temp[:, fcd_hist_temp[-1, :] ==
                                       emp_ks[-1, 0]]
        fcd_hist_mean = np.mean(fcd_histcum_nn, 1)
        ks_cost[t] = np.max(
            np.abs(fcd_hist_mean - emp_ks[:, 0]) / emp_ks[-1, 0])

    fcd_elapsed = time.time() - fcd_timestart
    print('Time using for cost function: ', fcd_elapsed)
    return ks_cost


def torch_corr(A):
    '''
    Self implemented correlation function used for GPU
    '''

    Amean = torch.mean(A, 1)
    Ax = A - torch.transpose(Amean.repeat(A.shape[1], 1), 0, 1)
    Astd = torch.mean(Ax**2, 1)
    Amm = torch.mm(Ax, torch.transpose(Ax, 0, 1)) / A.shape[1]
    Aout = torch.sqrt(torch.ger(Astd, Astd))
    Acor = Amm / Aout
    return Acor


def torch_corr2(A, B):
    '''
    Self implemented correlation function used for GPU
    '''

    Amean = torch.mean(A, 1)
    Ax = A - torch.transpose(Amean.repeat(A.shape[1], 1), 0, 1)
    Astd = torch.mean(Ax**2, 1)
    Bmean = torch.mean(B, 1)
    Bx = B - torch.transpose(Bmean.repeat(B.shape[1], 1), 0, 1)
    Bstd = torch.mean(Bx**2, 1)
    numerator = torch.mm(Ax, torch.transpose(Bx, 0, 1)) / A.shape[1]
    denominator = torch.sqrt(torch.ger(Astd, Bstd))
    torch_cor = numerator / denominator
    return torch_cor


def torch_arctanh(A):
    arctanh = 0.5 * torch.log((1 + A) / (1 - A))
    return arctanh


'''*************  Functions for computing FC & FCD costs   ****************'''


def CBIG_combined_cost_train(parameter, n_dup):
    '''
    This function is implemented to calcualted the FC correlation and FCD KS
    statistics combined cost for input parameter sets based on training data
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        n_dup:      number of times each parameter set is simulated
    Returns:
        total_cost: summation of FC correlation cost and FCD KS statistics cost
        corr_cost:  FC correlation cost
        ks_cost:    FCD KS statistics cost
    '''

    # Loading training data
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    emp_fcd = sio.loadmat('../../../input/Desikan_input/fcd_train.mat')
    emp_fcd = np.array(emp_fcd['train_aveM'])

    sc_mat_raw = csv_matrix_read('../../../input/Desikan_input/sc_train.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    emp_fc = csv_matrix_read('../../../input/Desikan_input/fc_train.csv')
    emp_fc = torch.from_numpy(emp_fc).type(torch.FloatTensor).cuda()

    # Calculating simualted BOLD signal using MFM
    bold_d = CBIG_mfm_multi_simulation(parameter, sc_mat, 14.4, n_dup)

    # Calculating FC correlation cost
    fc_cost = CBIG_FCcorrelation_multi_simulation(emp_fc, bold_d, n_dup)

    # Calculating FCD KS statistics cost
    fcd_cost = CBIG_FCDKSstat_multi_simulation(emp_fcd, bold_d, n_dup)

    # Calculating total cost
    total_cost = fc_cost + fcd_cost

    return total_cost, fc_cost, fcd_cost


def CBIG_combined_cost_validation(parameter, n_dup):
    '''
    This function is implemented to calcualted the FC correlation and FCD KS
    statistics combined cost for input parameter sets based on validation data
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        n_dup:      number of times each parameter set is simulated
    Returns:
        total_cost: summation of FC correlation cost and FCD KS statistics cost
        corr_cost:  FC correlation cost
        ks_cost:    FCD KS statistics cost
    '''

    # Loading validation data
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    emp_fcd = sio.loadmat('../../../input/Desikan_input/fcd_vali.mat')
    emp_fcd = np.array(emp_fcd['vali_aveM'])

    sc_mat_raw = csv_matrix_read('../../../input/Desikan_input/sc_vali.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    emp_fc = csv_matrix_read('../../../input/Desikan_input/fc_vali.csv')
    emp_fc = torch.from_numpy(emp_fc).type(torch.FloatTensor).cuda()

    # Calculating simualted BOLD signal using MFM
    bold_d = CBIG_mfm_multi_simulation(parameter, sc_mat, 14.4, n_dup)

    # Calculating FC correlation cost
    fc_cost = CBIG_FCcorrelation_multi_simulation(emp_fc, bold_d, n_dup)

    # Calculating FCD KS statistics cost
    fcd_cost = CBIG_FCDKSstat_multi_simulation(emp_fcd, bold_d, n_dup)

    # Calculating total cost
    total_cost = fc_cost + fcd_cost

    return total_cost, fc_cost, fcd_cost


def CBIG_combined_cost_test(parameter, n_dup):
    '''
    This function is implemented to calcualted the FC correlation and FCD KS
    statistics combined cost for input parameter sets based on test data
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        n_dup:      number of times each parameter set is simulated
    Returns:
        total_cost: summation of FC correlation cost and FCD KS statistics cost
        corr_cost:  FC correlation cost
        ks_cost:    FCD KS statistics cost
    '''

    # Loading validation data
    parameter = np.tile(parameter, [1, n_dup])
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    emp_fcd = sio.loadmat('../../../input/Desikan_input/fcd_test.mat')
    emp_fcd = np.array(emp_fcd['test_aveM'])

    sc_mat_raw = csv_matrix_read('../../../input/Desikan_input/sc_test.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    emp_fc = csv_matrix_read('../../../input/Desikan_input/fc_test.csv')
    emp_fc = torch.from_numpy(emp_fc).type(torch.FloatTensor).cuda()

    # Calculating simualted BOLD signal using MFM
    bold_d = CBIG_mfm_single_simulation(parameter, sc_mat, 14.4)

    # Calculating FC correlation cost
    fc_cost = CBIG_FCcorrelation_single_simulation(emp_fc, bold_d, n_dup)

    # Calculating FCD KS statistics cost
    fcd_cost = CBIG_FCDKSstat_single_simulation(emp_fcd, bold_d, n_dup)

    # Calculating total cost
    total_cost = fc_cost + fcd_cost

    return total_cost, fc_cost, fcd_cost


def CBIG_combined_cost_test_highres(parameter, n_dup):
    '''
    This function is implemented to calcualted the FC correlation and FCD KS
    statistics combined cost for input parameter sets based on test data
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        n_dup:      number of times each parameter set is simulated
    Returns:
        total_cost: summation of FC correlation cost and FCD KS statistics cost
        corr_cost:  FC correlation cost
        ks_cost:    FCD KS statistics cost
    '''

    # Loading validation data
    parameter = np.tile(parameter, [1, n_dup])
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    emp_fcd = sio.loadmat('../../../input/Desikan_input/fcd_test.mat')
    emp_fcd = np.array(emp_fcd['test_aveM'])

    sc_mat_raw = csv_matrix_read('../../../input/Desikan_input/sc_test.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    emp_fc = csv_matrix_read('../../../input/Desikan_input/fc_test.csv')
    emp_fc = torch.from_numpy(emp_fc).type(torch.FloatTensor).cuda()

    # Calculating simualted BOLD signal using MFM
    bold_d = CBIG_mfm_single_simulation(parameter, sc_mat, 14.4, d_t=0.001)

    # Calculating FC correlation cost
    fc_cost = CBIG_FCcorrelation_single_simulation(emp_fc, bold_d, n_dup)

    # Calculating FCD KS statistics cost
    fcd_cost = CBIG_FCDKSstat_single_simulation(emp_fcd, bold_d, n_dup)

    # Calculating total cost
    total_cost = fc_cost + fcd_cost

    return total_cost, fc_cost, fcd_cost


def CBIG_combined_cost_test_differwin(parameter, n_dup, window_indi):
    '''
    This function is implemented to calcualted the FC correlation and FCD KS
    statistics combined cost for input parameter sets based on test data
    Args:
        parameter:  (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength w
                    parameter[N:2*N]: external input I
                    parameter[2*N]: Gloable constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        n_dup:      number of times each parameter set is simulated
        window_indi:determine the size the sliding window size
                    'low':  window size is 43
                    'high': window size is 125
    Returns:
        total_cost: summation of FC correlation cost and FCD KS statistics cost
        corr_cost:  FC correlation cost
        ks_cost:    FCD KS statistics cost
    '''

    # Loading validation data
    parameter = np.tile(parameter, [1, n_dup])
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    if window_indi == 'low':
        emp_fcd = sio.loadmat(
            '../../../input/Desikan_input/fcd_test_low_window.mat')
        window = 43
    elif window_indi == 'high':
        emp_fcd = sio.loadmat(
            '../../../input/Desikan_input/fcd_test_high_window.mat')
        window = 125
    else:
        raise ValueError('Input is not acceptable.')
    emp_fcd = np.array(emp_fcd['test_aveM'])

    sc_mat_raw = csv_matrix_read('../../../input/Desikan_input/sc_test.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    emp_fc = csv_matrix_read('../../../input/Desikan_input/fc_test.csv')
    emp_fc = torch.from_numpy(emp_fc).type(torch.FloatTensor).cuda()

    # Calculating simualted BOLD signal using MFM
    bold_d = CBIG_mfm_single_simulation(parameter, sc_mat, 14.4)

    # Calculating FC correlation cost
    fc_cost = CBIG_FCcorrelation_single_simulation(emp_fc, bold_d, n_dup)

    # Calculating FCD KS statistics cost
    fcd_cost = CBIG_FCDKSstat_single_simulation(
        emp_fcd, bold_d, n_dup, window_size=window)

    # Calculating total cost
    total_cost = fc_cost + fcd_cost

    return total_cost, fc_cost, fcd_cost


def csv_matrix_read(filename):

    csv_file = open(filename, "r")
    read_handle = csv.reader(csv_file)
    out_list = []
    R = 0
    for row in read_handle:
        out_list.append([])
        for col in row:
            out_list[R].append(float(col))
        R = R + 1
    out_array = np.array(out_list)
    csv_file.close()
    return out_array
