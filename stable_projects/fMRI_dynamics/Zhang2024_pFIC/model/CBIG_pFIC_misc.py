# /usr/bin/env python
"""
Written by Shaoshi Zhang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import csv
import math
import time
import warnings
import torch
import numpy as np
import scipy.io as sio
from scipy.optimize import fsolve


def CBIG_mfm_multi_simulation(parameter,
                              sc_mat,
                              simulation_period,
                              n_dup,
                              d_t,
                              low_mem,
                              parser,
                              rE_lower_bound=2.7,
                              rE_upper_bound=3.3):
    """

    Generate simulated BOLD signals with different candidate parameters.
    This function is mainly used during the model training stage.
    Each set (column) of the candidate parameters generates multiple
    simulated BOLD signals (n_dup times, to be exact), each time with a
    different noise instantiation. But the corresponding repetition of
    different candidate parameters share the same noise instantiation
    (in general, noises of the x-th simulation of the y-th set of
    candidate parameters and noises of the x-th simulation of the
    z-th set of candidate parameters are the same).
    Input:
        -parameter: #feature x #set candidate parameter matrix.
                    Each row corresponds to a feature (for example,
                    wEE of the n-th ROI or global constant G)
        -sc_mat: #ROI x #ROI matrix, structural connectivity
        -simulation_length: simulation period in minutes
                    (for example, 14.4 for HCP young adult dataset)
        -n_dup: number of repetition simulated for each set of
                    candidate parameters
        -d_t: temporal step size for simulation in seconds
                    (for example, 0.006)
        -low_mem: an option to save memory when d_t is too small
                    (generate noise on the fly)
        -parser: a parsed configuration file specifying model parameters
                    (see example.ini)
        -rE_lower_bound, rE_upper_bound: lower and upper bounds of
                    allowed excitatory firing rate thresholds,
                    simulated BOLD with firing rate higher or
                    lower than the thresholds will be set
                    to NaN and excluded
    Output:
        -bold: #ROI x (#set x #repetition) x #TR 3D matrix.
               Simulated BOLD signal. The second dimension is the product
               of the number of candidate parameter sets and number of
               repeated simulations (i.e. n_dup)
        -S_E_all: #ROI x (#set x  #repetition) 2D matrix. Simulated
                excitatory synaptic gating variables. Temporally averaged.
        -S_I_all: #ROI x (#set x #repetition) 2D matrix. Simulated
                inhibitory synaptic gating variables. Temporally averaged.
        -r_E_all: #ROI x (#set x #repetition) 2D matrix. Simulated excitatory
                firing rates. Temporally averaged.
        -J_I: #ROI x (#set x #repetition) 2D matrix.
                Feedback inhibition strength. Equivalent to wIE.

    """
    torch.set_default_tensor_type('torch.cuda.DoubleTensor')

    # Initializing system parameters
    warmup = int(
        parser['system']['warmup'])  # number of frames used for warmup
    kstart = 0.
    # pre-simulation period in minutes, this part will be excluded
    # (together with warmup) in the end
    t_pre = 60 * float(parser['system']['t_pre'])
    kend = t_pre + 60 * simulation_period

    # Setting sampling ratio
    k_p = torch.arange(kstart, kend + d_t, d_t)  # number of simulated frames
    n_set = parameter.shape[1]  # number of sets of candidate parameters
    n_num = n_dup * n_set  # total number of simulations
    parameter = parameter.repeat(1, n_dup)
    n_node = sc_mat.shape[0]  # number of ROIs
    n_frame = k_p.shape[0]
    TR = float(parser['BOLD']['TR'])

    # Initializing neural and hemodynamic activity
    S_E = torch.zeros((n_node, n_num))
    S_I = torch.zeros((n_node, n_num))
    I_I_ave = float(parser['inhibitory']['I_I_ss']) * torch.ones(n_node, n_num)
    f_mat = torch.ones((n_node, n_num, 4))
    z_t = torch.zeros((n_node, n_num))
    f_t = torch.ones((n_node, n_num))
    v_t = torch.ones((n_node, n_num))
    q_t = torch.ones((n_node, n_num))
    f_mat[:, :, 0] = z_t

    S_E[:, :] = float(parser['excitatory']['S_E_ss'])  # initialization
    S_I[:, :] = float(parser['inhibitory']['S_I_ss'])  # initialization
    sigma = parameter[2 * n_node + 1:3 * n_node + 1, :]  # noise amplitude
    w_l = k_p.shape[0]
    if low_mem == 0:
        d_w = math.sqrt(d_t) * torch.randn(n_dup, n_node,
                                           w_l + warmup)  # make some noise
    p_costant = float(parser['hemodynamic']['p_constant'])
    v_0 = float(parser['hemodynamic']['v_0'])
    k_1 = 4.3 * 28.265 * 3 * 0.0331 * p_costant
    k_2 = 0.47 * 110 * 0.0331 * p_costant
    k_3 = 0.53
    count = 0
    y_bold = torch.zeros((n_node, n_num, int(n_frame / (TR / d_t) + 1)))

    # Parameters for firing rate
    a_I = float(parser['inhibitory']['a_I'])
    b_I = float(parser['inhibitory']['b_I'])
    d_I = float(parser['inhibitory']['d_I'])

    # solve I_I_ave
    S_E_ave = float(parser['excitatory']['S_E_ss'])
    I_E_ave = float(parser['excitatory']['I_E_ss'])
    J_NMDA = float(parser['excitatory']['J_NMDA'])
    w_EE = parameter[0:n_node, :]
    G = parameter[2 * n_node, :]
    w_EI = parameter[n_node:2 * n_node, :]
    W_E = float(parser['excitatory']['W_E'])
    I0 = float(parser['neuralMassModel']['I0'])
    tau_I = float(parser['inhibitory']['tau_I'])
    w_EI = w_EI.cpu().numpy()
    I_I_ave = I_I_ave.cpu().numpy()
    # Analytically compute feedback inhibition (wIE). See Demirtas 2019 Neuron
    # for more details
    for i in range(n_num):
        I_I_ave_one_set = np.atleast_2d(I_I_ave[:, i]).T
        w_EI_one_set = np.atleast_2d(w_EI[:, i]).T
        I_I_ave[:, i], infodict, ier, mseg = fsolve(
            I_I_fixed_pt,
            I_I_ave_one_set,
            args=(w_EI_one_set, parser),
            full_output=True)

    I_I_ave = torch.from_numpy(I_I_ave).type(torch.DoubleTensor).cuda()
    S_I_ave = tau_I * (a_I * I_I_ave - b_I) / (
        1 - torch.exp(-d_I * (a_I * I_I_ave - b_I)))

    # calculate J_I (or wIE)
    J_I = torch.div(
        W_E * I0 + w_EE * J_NMDA * S_E_ave + G * J_NMDA * torch.sum(
            sc_mat, 1).view(-1, 1).repeat(1, n_num) * S_E_ave - I_E_ave,
        S_I_ave)

    # Warm up
    start = time.time()
    print('Warm-up simulation starts...')
    r_E_all = torch.zeros(n_node, n_num)
    S_E_all = torch.zeros(n_node, n_num)
    S_I_all = torch.zeros(n_node, n_num)
    for i in range(warmup):
        dS_E, dS_I, r_E = CBIG_mfm_rfMRI_ode(S_E, S_I, J_I, parameter, sc_mat,
                                             parser)
        if low_mem == 0:
            noise = d_w[:, :, i].repeat(1, 1, n_set).contiguous().view(
                -1, n_node)
            S_E = S_E + dS_E * d_t + sigma * torch.transpose(noise, 0, 1)
            S_I = S_I + dS_I * d_t + sigma * torch.transpose(noise, 0, 1)
        else:
            S_E = S_E + dS_E * d_t + sigma * torch.randn(n_node, n_dup).repeat(
                1, n_set) * math.sqrt(d_t)
            S_I = S_I + dS_I * d_t + sigma * torch.randn(n_node, n_dup).repeat(
                1, n_set) * math.sqrt(d_t)

    print('Simulation starts...')
    for i in range(n_frame):
        dS_E, dS_I, r_E = CBIG_mfm_rfMRI_ode(S_E, S_I, J_I, parameter, sc_mat,
                                             parser)
        if low_mem == 0:
            noise = d_w[:, :, i].repeat(1, 1, n_set).contiguous().view(
                -1, n_node)
            S_E = S_E + dS_E * d_t + sigma * torch.transpose(noise, 0, 1)
            S_I = S_I + dS_I * d_t + sigma * torch.transpose(noise, 0, 1)
        else:
            S_E = S_E + dS_E * d_t + sigma * torch.randn(n_node, n_dup).repeat(
                1, n_set) * math.sqrt(d_t)
            S_I = S_I + dS_I * d_t + sigma * torch.randn(n_node, n_dup).repeat(
                1, n_set) * math.sqrt(d_t)

        d_f = CBIG_mfm_rfMRI_BW_ode(S_E, f_mat, parser)
        f_mat = f_mat + d_f * d_t
        z_t, f_t, v_t, q_t = torch.chunk(f_mat, 4, dim=2)
        y_bold_temp = 100 / p_costant * v_0 * (
            k_1 * (1 - q_t) + k_2 * (1 - q_t / v_t) + k_3 * (1 - v_t))
        # downsample to match TR of BOLD signals
        y_bold[:, :, count] = y_bold_temp[:, :, 0]
        count = count + ((i + 1) % (int(round(TR / d_t))) == 0) * 1

        r_E_all = r_E_all + r_E
        S_E_all = S_E_all + S_E
        S_I_all = S_I_all + S_I

    elapsed = time.time() - start
    r_E_all = r_E_all / n_frame
    S_E_all = S_E_all / n_frame
    S_I_all = S_I_all / n_frame

    # if the excitatory firing rate of any node is below rE_lower_bound or
    # above rE_upper_bound, the first node of the last frame of BOLD is set
    # to NaN (essentially render this set of parameters un-usable)
    for i in range(n_num):
        if (rE_lower_bound > r_E_all[:, i]).any() or (r_E_all[:, i] >
                                                      rE_upper_bound).any():
            y_bold[:, i, :] = float('nan')

    # Remove pre-simulation period
    cut_index = int(t_pre / TR)
    bold = y_bold[:, :, cut_index + 1:y_bold.shape[2]]

    print('End of simulation. Time used:', elapsed, 'seconds')
    return bold, S_E_all, S_I_all, r_E_all, J_I


def CBIG_mfm_single_simulation(parameter,
                               sc_mat,
                               simulation_period,
                               d_t,
                               parser,
                               rE_lower_bound=2.7,
                               rE_upper_bound=3.3):
    """

    Generate simulated BOLD signals with different candidate parameters.
    This function is mainly used during the model testing stage. Each set (
    column) of the candidate parameters generates only 1 set of simulated
    BOLD signals.

    Input:
        -parameter: #feature x #set candidate parameter matrix. Each row
            corresponds to a feature (for example, wEE of the n-th
            ROI or global constant G)
        -sc_mat: #ROI x #ROI matrix, structural connectivity
        -simulation_length: simulation period in minutes
            (for example, 14.4 for HCP young adult dataset)
        -d_t: temporal step size for simulation in seconds,
            normally a shorter step size than during training
            (for example, 0.0005)
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
        -rE_lower_bound, rE_upper_bound: lower and upper bounds of allowed
            excitatory firing rate thresholds,
            simulated BOLD with firing rate higher or lower than the
            thresholds will be set to NaN and excluded
    Output:
        -bold: #ROI x (#set x #repetition) x #TR 3D matrix.
            Simulated BOLD signal. The second dimension is the product
            of the number of candidate parameter sets and number of
            repeated simulations (i.e. n_dup)
        -S_E_all: #ROI x (#set x  #repetition) 2D matrix. Simulated
            excitatory synaptic gating variables. Temporally averaged.
        -S_I_all: #ROI x (#set x #repetition) 2D matrix. Simulated
            inhibitory synaptic gating variables. Temporally averaged.
        -r_E_all: #ROI x (#set x #repetition) 2D matrix. Simulated
            excitatory firing rates. Temporally averaged.
        -J_I: #ROI x (#set x #repetition) 2D matrix. Feedback
            inhibition strength. Equivalent to wIE.

    """

    torch.set_default_tensor_type('torch.cuda.DoubleTensor')

    # Initializing system parameters
    warmup = int(
        parser['system']['warmup'])  # number of frames used for warmup
    kstart = 0.
    # pre-simulation period in minutes, this part will be excluded
    # (together with warmup) in the end
    t_pre = 60 * float(parser['system']['t_pre'])
    kend = t_pre + 60 * simulation_period
    TR = float(parser['BOLD']['TR'])

    # sampling ratio
    k_p = torch.arange(kstart, kend + d_t, d_t)
    n_node = sc_mat.shape[0]
    n_frame = k_p.shape[0]
    n_num = parameter.shape[1]

    # # Initializing neural and hemodynamic activity
    S_E = torch.zeros((n_node, n_num))
    S_I = torch.zeros((n_node, n_num))
    I_I_ave = float(parser['inhibitory']['I_I_ss']) * torch.ones(n_node, n_num)
    f_mat = torch.ones((n_node, n_num, 4))
    z_t = torch.zeros((n_node, n_num))
    f_t = torch.ones((n_node, n_num))
    v_t = torch.ones((n_node, n_num))
    q_t = torch.ones((n_node, n_num))
    f_mat[:, :, 0] = z_t

    S_E[:, :] = float(parser['excitatory']['S_E_ss'])  # initialization
    S_I[:, :] = float(parser['inhibitory']['S_I_ss'])  # initialization
    sigma = parameter[2 * n_node + 1:3 * n_node + 1, :]
    p_costant = float(parser['hemodynamic']['p_constant'])
    v_0 = float(parser['hemodynamic']['v_0'])
    k_1 = 4.3 * 28.265 * 3 * 0.0331 * p_costant
    k_2 = 0.47 * 110 * 0.0331 * p_costant
    k_3 = 0.53
    count = 0
    y_bold = torch.zeros((n_node, n_num, int(n_frame / (TR / d_t) + 1)))

    # Parameters for firing rate
    a_I = float(parser['inhibitory']['a_I'])
    b_I = float(parser['inhibitory']['b_I'])
    d_I = float(parser['inhibitory']['d_I'])

    # solve I_I_ave
    S_E_ave = float(parser['excitatory']['S_E_ss'])
    I_E_ave = float(parser['excitatory']['I_E_ss'])
    J_NMDA = float(parser['excitatory']['J_NMDA'])
    w_EE = parameter[0:n_node, :]
    G = parameter[2 * n_node, :]
    w_EI = parameter[n_node:2 * n_node, :]
    W_E = float(parser['excitatory']['W_E'])
    I0 = float(parser['neuralMassModel']['I0'])
    tau_I = float(parser['inhibitory']['tau_I'])
    w_EI = w_EI.cpu().numpy()
    I_I_ave = I_I_ave.cpu().numpy()
    # See Demirtas 2019 for more details
    for i in range(n_num):
        I_I_ave_one_set = np.atleast_2d(I_I_ave[:, i]).T
        w_EI_one_set = np.atleast_2d(w_EI[:, i]).T
        I_I_ave[:, i], infodict, ier, mseg = fsolve(
            I_I_fixed_pt,
            I_I_ave_one_set,
            args=(w_EI_one_set, parser),
            full_output=True)
    I_I_ave = torch.from_numpy(I_I_ave).type(torch.DoubleTensor).cuda()
    S_I_ave = tau_I * (a_I * I_I_ave - b_I) / (
        1 - torch.exp(-d_I * (a_I * I_I_ave - b_I)))

    # calculate J_I (or wIE)
    J_I = torch.div(
        W_E * I0 + w_EE * J_NMDA * S_E_ave + G * J_NMDA * torch.sum(
            sc_mat, 1).view(-1, 1).repeat(1, n_num) * S_E_ave - I_E_ave,
        S_I_ave)

    # Warm up
    start = time.time()
    print('Warm-up simulation starts...')
    r_E_all = torch.zeros(n_node, n_num)
    S_E_all = torch.zeros(n_node, n_num)
    S_I_all = torch.zeros(n_node, n_num)
    for i in range(warmup):
        dS_E, dS_I, r_E = CBIG_mfm_rfMRI_ode(S_E, S_I, J_I, parameter, sc_mat,
                                             parser)
        S_E = S_E + dS_E * d_t + sigma * torch.randn(n_node,
                                                     n_num) * math.sqrt(d_t)
        S_I = S_I + dS_I * d_t + sigma * torch.randn(n_node,
                                                     n_num) * math.sqrt(d_t)

    # Main body
    print('Simulation starts...')
    for i in range(n_frame):
        dS_E, dS_I, r_E = CBIG_mfm_rfMRI_ode(S_E, S_I, J_I, parameter, sc_mat,
                                             parser)
        S_E = S_E + dS_E * d_t + sigma * torch.randn(n_node,
                                                     n_num) * math.sqrt(d_t)
        S_I = S_I + dS_I * d_t + sigma * torch.randn(n_node,
                                                     n_num) * math.sqrt(d_t)

        d_f = CBIG_mfm_rfMRI_BW_ode(S_E, f_mat, parser)
        f_mat = f_mat + d_f * d_t
        z_t, f_t, v_t, q_t = torch.chunk(f_mat, 4, dim=2)
        y_bold_temp = 100 / p_costant * v_0 * (
            k_1 * (1 - q_t) + k_2 * (1 - q_t / v_t) + k_3 * (1 - v_t))
        y_bold[:, :, count] = y_bold_temp[:, :, 0]
        count = count + ((i + 1) % (int(round(TR / d_t))) == 0) * 1

        r_E_all = r_E_all + r_E
        S_E_all = S_E_all + S_E
        S_I_all = S_I_all + S_I

    elapsed = time.time() - start
    r_E_all = r_E_all / n_frame
    S_E_all = S_E_all / n_frame
    S_I_all = S_I_all / n_frame

    # if the excitatory firing rate of any node is below rE_lower_bound
    # or above rE_upper_bound, the first node of the last frame of BOLD is
    # set to NaN (essentially render this set of parameters un-usable)
    for i in range(n_num):
        if (rE_lower_bound > r_E_all[:, i]).any() or (r_E_all[:, i] >
                                                      rE_upper_bound).any():
            y_bold[:, i, :] = float('nan')

    # Remove pre-simulation period
    cut_index = int(t_pre / TR)
    bold = y_bold[:, :, cut_index + 1:y_bold.shape[2]]

    print('End of simulation. Time used:', elapsed, 'seconds')

    return bold, S_E_all, S_I_all, r_E_all, J_I


def CBIG_mfm_rfMRI_ode(S_E, S_I, J_I, parameter, sc_mat, parser):
    """

    The function implements Deco et al., 2014 ODE of
    feedback inhibition control model
    Input:
        -S_E: excitatory synaptic gating variable of the preceding frame
        -S_I: inhibitory synaptic gating variable of the preceding frame
        -J_I: feedback inhibition control strength (or wIE)
        -parameter: model parameter (wEE, wEI, G, sigma)
        -sc_mat: structural connectivity matrix
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
    Output:
        -dS_E: temporal derivative of excitatory synaptic gating variable
            of the current frame
        -dS_I: temporal derivative of inhibitory synaptic gating variable
            of the current frame
        -r_E: excitatory firing rate of the current frame
    """

    torch.set_default_tensor_type('torch.cuda.DoubleTensor')

    n_node = sc_mat.shape[0]

    # current parameters
    J_NMDA = float(parser['excitatory']['J_NMDA'])
    w_EE = parameter[0:n_node, :]
    w_EI = parameter[n_node:2 * n_node, :]
    G = parameter[2 * n_node, :]
    W_E = float(parser['excitatory']['W_E'])
    W_I = float(parser['inhibitory']['W_I'])
    I0 = float(parser['neuralMassModel']['I0'])

    # firing rate parameters
    a_E = float(parser['excitatory']['a_E'])
    b_E = float(parser['excitatory']['b_E'])
    d_E = float(parser['excitatory']['d_E'])

    a_I = float(parser['inhibitory']['a_I'])
    b_I = float(parser['inhibitory']['b_I'])
    d_I = float(parser['inhibitory']['d_I'])

    # synaptic gating variable parameters
    tau_E = float(parser['excitatory']['tau_E'])
    tau_I = float(parser['inhibitory']['tau_I'])
    gamma = float(parser['neuralMassModel']['gamma'])

    # current
    I_E = W_E * I0 + J_NMDA * w_EE * S_E + J_NMDA * G.repeat(
        n_node, 1) * torch.mm(sc_mat, S_E) - J_I * S_I
    I_I = W_I * I0 + J_NMDA * w_EI * S_E - S_I

    # firing rate
    r_E = (a_E * I_E - b_E) / (1 - torch.exp(-d_E * (a_E * I_E - b_E)))
    r_I = (a_I * I_I - b_I) / (1 - torch.exp(-d_I * (a_I * I_I - b_I)))

    # temporal derivative of synaptic gating variable (noise is added at the
    # forward simulation step)
    dS_E = -S_E / tau_E + (1 - S_E) * gamma * r_E
    dS_I = -S_I / tau_I + r_I

    return dS_E, dS_I, r_E


def I_I_fixed_pt(I_I_ave, w_EI, parser):
    """

    This function calculates the steady-state inhibitory current based on the
      given wEI. See Demirtas et al., 2019 for more details. This function is
      integrated inside a fsolve function for an analytic solution of wIE
    Input:
        -I_I_ave: an initial guess of the steady-state inhibitory current
        -w_EI: exictatory to inhibitory strength
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
    Output:
        -I_I_ave_fixed_point.T.flatten(): solved steady-state inhibitory
            current

    """

    I_I_ave = np.atleast_2d(I_I_ave).T
    W_I = float(parser['inhibitory']['W_I'])
    J_NMDA = float(parser['excitatory']['J_NMDA'])
    S_E_ave = float(parser['excitatory']['S_E_ss'])
    a_I = float(parser['inhibitory']['a_I'])
    b_I = float(parser['inhibitory']['b_I'])
    d_I = float(parser['inhibitory']['d_I'])
    tau_I = float(parser['inhibitory']['tau_I'])
    I0 = float(parser['neuralMassModel']['I0'])

    I_I_ave_fixed_point = -I_I_ave + W_I * I0 + J_NMDA * w_EI * S_E_ave - (
        a_I * I_I_ave - b_I) / (1 - np.exp(-d_I *
                                           (a_I * I_I_ave - b_I))) * tau_I
    return I_I_ave_fixed_point.T.flatten()


def CBIG_mfm_rfMRI_BW_ode(y_t, F, parser):
    """

    This function implements the Balloon-Windekessel model that converts
    neuronal activities to BOLD signals. See Wang 2019 (SciAdv) and Stephan
    2007 for more details
    Input:
        -y_t: neuronal activities. This should be S_E for this model.
        -F: #ROI x #sets of parameters x 4 matrix. Each matrix along the
            3rd dimension represents a hemodynamic variable
            in this order: z, f, v, q
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
    Output:
        -dF: temporal derivatives of F

    """

    torch.set_default_tensor_type('torch.cuda.DoubleTensor')

    # Hemodynamic model parameters
    beta = float(parser['hemodynamic']['beta'])
    gamma = float(parser['hemodynamic']['gamma'])
    tau = float(parser['hemodynamic']['tau'])
    alpha = float(parser['hemodynamic']['alpha'])
    p = float(parser['hemodynamic']['p_constant'])
    n_nodes = y_t.shape[0]
    n_set = y_t.shape[1]
    z = F[:, :, 0]
    f = F[:, :, 1]
    v = F[:, :, 2]
    q = F[:, :, 3]

    # Calculate derivatives
    dF = torch.zeros((n_nodes, n_set, 4))
    dF[:, :, 0] = y_t - beta * z - gamma * (f - 1)
    dF[:, :, 1] = z
    dF[:, :, 2] = 1 / tau * (f - pow(v, (1 / alpha)))
    dF[:, :, 3] = 1 / tau * (f / p * (1 - pow(
        (1 - p), (1 / f))) - q / v * pow(v, (1 / alpha)))
    return dF


def CBIG_FC_multi_simulation(emp_fc, bold, n_dup, FC_weights):
    """

    This function calculates FC correlation cost and L1 cost between
    simulated and empirical FC
    Input:
        -emp_fc: #ROI x #ROI empirical FC matrix
        -bold: #ROI x (#set x #repetition) x #TR 3D matrix.
            Simulated BOLD signal. Output from CBIG_mfm_multi_simulation
        -n_dup: number of repetition for each set of candidate parameters
        -FC_weights: 2 x 1 array. The first value is FC correlation loss
             weight, the second value is FC L1 loss weight
    Output:
        -FC_corr_weight * corr_cost: weighted FC correlation loss
        -FC_L1_weight * L1_cost: weighted FC L1 loss

    """

    fc_timestart = time.time()
    FC_corr_weight = FC_weights[0]
    FC_L1_weight = FC_weights[1]

    # vectorize simulated FC
    n_num = bold.shape[1]
    n_set = int(n_num / n_dup)
    n_node = emp_fc.shape[0]
    fc_mask = torch.triu(torch.ones(n_node, n_node), 1) == 1
    vect_len = int(n_node * (n_node - 1) / 2)
    sim_fc_vector = torch.zeros(n_num, vect_len)
    for i in range(n_num):
        sim_fc = torch_corr(bold[:, i, :])
        sim_fc_vector[i, :] = sim_fc[fc_mask]

    # Average the simulated FCs within the same parameter set

    # (across different noise instantiations)
    sim_fc_vector[sim_fc_vector != sim_fc_vector] = 0
    sim_fc_numerator = torch.zeros(n_set, vect_len)
    sim_fc_denom = torch.zeros(n_set, 1)
    for k in range(n_dup):
        sim_fc_numerator = sim_fc_numerator + sim_fc_vector[k * n_set:
                                                            (k + 1) * n_set, :]
        sim_fc_denom = sim_fc_denom + (
            sim_fc_vector[k * n_set:(k + 1) * n_set, 0:1] != 0).float()

    sim_fc_denom[sim_fc_denom == 0] = np.nan
    sim_fc_ave = sim_fc_numerator / sim_fc_denom

    # FC correlation cost
    emp_fcm = emp_fc[fc_mask].repeat(n_set, 1)
    corr_mass = torch_corr2(sim_fc_ave, emp_fcm)
    corr_cost = torch.diag(corr_mass)
    corr_cost = corr_cost.cpu().numpy()
    corr_cost = 1 - corr_cost
    corr_cost[np.isnan(
        corr_cost
    )] = 5  # an arbitrarily high value to indicate NaN FC correlation cost

    # FC L1 cost
    L1_cost = torch.abs(torch.mean(sim_fc_ave, 1) - torch.mean(emp_fcm, 1))
    L1_cost = L1_cost.cpu().numpy()
    L1_cost[np.isnan(
        L1_cost)] = 5  # an arbitrarily high value to indicate NaN FC L1 cost

    fc_elapsed = time.time() - fc_timestart
    print('Time taken to calculate FC cost:', fc_elapsed, 'seconds')

    return FC_corr_weight * corr_cost, FC_L1_weight * L1_cost


def CBIG_FC_single_simulation(emp_fc, bold, n_dup, FC_weights):
    """

    This function calculates FC correlation cost and L1 cost between
    simulated and empirical FC
    Input:
        -emp_fc: #ROI x #ROI empirical FC matrix
        -bold: #ROI x (#set x #repetition) x #TR 3D matrix.
            Simulated BOLD signal. Output from CBIG_mfm_single_simulation
        -FC_weights: 2 x 1 array. The first value is FC correlation loss
            weight, the second value is FC L1 loss weight
    Output:
        -FC_corr_weight * corr_cost: weighted FC correlation loss
        -FC_L1_weight * L1_cost: weighted FC L1 loss

    """

    fc_timestart = time.time()
    FC_corr_weight = FC_weights[0]
    FC_L1_weight = FC_weights[1]

    # Vectorize simulated FC
    n_num = bold.shape[1]
    n_node = emp_fc.shape[0]
    fc_mask = torch.triu(torch.ones(n_node, n_node), 1) == 1
    vect_len = int(n_node * (n_node - 1) / 2)
    sim_fc_vector = torch.zeros(n_num, vect_len)
    for i in range(n_num):
        sim_fc = torch_corr(bold[:, i, :])
        sim_fc_vector[i, :] = sim_fc[fc_mask]

    sim_fc_numpy = sim_fc_vector.cpu().numpy()
    emp_fc_numpy = emp_fc[fc_mask].cpu().numpy()
    time_dup = int(n_num / n_dup)
    corr_cost = np.zeros(time_dup)
    L1_cost = np.zeros(time_dup)

    for t in range(time_dup):
        sim_fc_numpy_temp = sim_fc_numpy[t * n_dup:(t + 1) * n_dup, :]
        sim_fc_mean = np.nanmean(sim_fc_numpy_temp, 0)
        corrmean_temp = np.corrcoef(sim_fc_mean, emp_fc_numpy)
        corr_cost[t] = 1 - corrmean_temp[1, 0]  # FC correlation loss
        L1_cost[t] = np.abs(
            np.mean(sim_fc_mean) - np.mean(emp_fc_numpy))  # FC L1 loss

    fc_elapsed = time.time() - fc_timestart
    print('Time taken to calculate FC cost:', fc_elapsed, 'seconds')
    return FC_corr_weight * corr_cost, FC_L1_weight * L1_cost


def CBIG_FCD_multi_simulation(emp_fcd_cdf, bold, n_dup, FCD_weight, parser):
    """

    This function calculates FCD KS statistics between
    simulated and empirical FCD CDF
    Input:
        -emp_fcd_cdf: 1 x 10000 empirical FCD CDF
        -bold: #ROI x (#set x #repetition) x #TR 3D matrix.
            Simulated BOLD signal. Output from CBIG_mfm_multi_simulation
        -n_dup: number of repetition for each set of candidate parameters
        -FCD_weights: 1 x 1 scalar
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
    Output:
        -FCD_weight * ks_cost: weighted FCD KS loss

    """

    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    fcd_timestart = time.time()

    # Initializing the FC and FCD masks
    n_num = bold.shape[1]
    n_set = int(n_num / n_dup)
    n_node = bold.shape[0]
    window_size = int(parser['system']['window_size'])
    n_frame = int(parser['BOLD']['n_frame'])
    time_lengh = n_frame - window_size + 1
    sub_num = 10
    resid_num = n_num % sub_num
    fc_edgenum = int(n_node * (n_node - 1) / 2)
    fc_mask = torch.triu(torch.ones(n_node, n_node), 1) == 1
    fc_maskm = torch.zeros(n_node * sub_num,
                           n_node * sub_num).type(torch.cuda.ByteTensor)
    for i in range(sub_num):
        fc_maskm[n_node * i:n_node * (i + 1), n_node * i:n_node * (i + 1)] = \
            fc_mask

    fc_mask_resid = torch.zeros(n_node * resid_num,
                                n_node * resid_num).type(torch.cuda.ByteTensor)
    for i in range(resid_num):
        j = i + 1
        fc_mask_resid[n_node * i:n_node * j, n_node * i:n_node * j] = fc_mask

    fcd_mask = torch.triu(torch.ones(time_lengh, time_lengh), 1) == 1

    # Calculate simulated FCD CDF
    fcd_hist = np.ones([10000, n_num])
    fc_mat = torch.zeros(fc_edgenum, sub_num, time_lengh)
    batch_num = math.floor(n_num / sub_num)
    fc_resid = torch.zeros(fc_edgenum, resid_num, time_lengh)

    for b in range(batch_num):
        bold_temp = bold[:, b * sub_num:(b + 1) * sub_num, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, n_frame)
        for i in range(0, time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_maskm]
            fc_mat[:, :, i] = torch.transpose(
                cor_temp.view(sub_num, fc_edgenum), 0, 1)

        for j in range(0, sub_num):
            fcd_temp = torch_corr(torch.transpose(fc_mat[:, j, :], 0, 1))
            fcd_hist_temp = np.histogram(
                fcd_temp[fcd_mask].cpu().numpy(), bins=10000, range=(-1., 1.))
            fcd_hist[:, j + b * sub_num] = fcd_hist_temp[0]

    if resid_num != 0:
        bold_temp = bold[:, batch_num * sub_num:n_num, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, n_frame)
        for i in range(time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_mask_resid]
            fc_resid[:, :, i] = torch.transpose(
                cor_temp.view(resid_num, fc_edgenum), 0, 1)

        for j in range(resid_num):
            fcd_temp = torch_corr(torch.transpose(fc_resid[:, j, :], 0, 1))
            fcd_hist_temp = np.histogram(
                fcd_temp[fcd_mask].cpu().numpy(), bins=10000, range=(-1., 1.))
            fcd_hist[:, j + sub_num * batch_num] = fcd_hist_temp[0]

    fcd_histcum = np.cumsum(fcd_hist, 0)
    fcd_histcumM = fcd_histcum.copy()
    fcd_histcumM[:, fcd_histcum[-1, :] != emp_fcd_cdf[-1, 0]] = 0

    # Calculate KS  cost
    fcd_histcum_temp = np.zeros((10000, n_set))
    fcd_histcum_num = np.zeros((1, n_set))
    for k in range(n_dup):
        fcd_histcum_temp = fcd_histcum_temp + fcd_histcumM[:, k * n_set:
                                                           (k + 1) * n_set]
        fcd_histcum_num = fcd_histcum_num + (
            fcd_histcumM[-1, k * n_set:(k + 1) * n_set] == emp_fcd_cdf[-1, 0])
    fcd_histcum_ave = fcd_histcum_temp / fcd_histcum_num
    ks_diff = np.abs(fcd_histcum_ave - np.tile(emp_fcd_cdf, [1, n_set]))
    ks_cost = ks_diff.max(0) / emp_fcd_cdf[-1, 0]
    ks_cost[fcd_histcum_ave[-1, :] != emp_fcd_cdf[-1, 0]] = 5

    fcd_elapsed = time.time() - fcd_timestart
    print('Time taken to calculate FCD cost:', fcd_elapsed, 'seconds')
    return FCD_weight * ks_cost


def CBIG_FCD_single_simulation(emp_ks, bold_d, n_dup, FCD_weight, parser):
    """

    This function calculates FCD KS statistics between
    simulated and empirical FCD CDF
    Input:
        -emp_fcd_cdf: 1 x 10000 empirical FCD CDF
        -bold: #ROI x (#set x #repetition) x #TR 3D matrix.
            Simulated BOLD signal. Output from CBIG_mfm_single_simulation
        -n_dup: number of repetition for each set of candidate parameters
        -FCD_weights: 1 x 1 scalar
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
    Output:
        -FCD_weight * ks_cost: weighted FCD KS loss

    """

    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    fcd_timestart = time.time()

    # Initializing the FC and FCD masks
    n_set = bold_d.shape[1]
    n_node = bold_d.shape[0]
    window_size = int(parser['system']['window_size'])
    n_frame = int(parser['BOLD']['n_frame'])
    time_lengh = n_frame - window_size + 1
    sub_num = 10
    resid_num = n_set % sub_num
    fc_edgenum = int(n_node * (n_node - 1) / 2)
    fc_mask = torch.triu(torch.ones(n_node, n_node), 1) == 1
    fc_maskm = torch.zeros(n_node * sub_num,
                           n_node * sub_num).type(torch.cuda.ByteTensor)

    for i in range(sub_num):
        fc_maskm[n_node * i:n_node * (i + 1), n_node * i:n_node *
                 (i + 1)] = fc_mask

    fc_mask_resid = torch.zeros(n_node * resid_num,
                                n_node * resid_num).type(torch.cuda.ByteTensor)
    for i in range(resid_num):
        j = i + 1
        fc_mask_resid[n_node * i:n_node * j, n_node * i:n_node * j] = fc_mask

    fcd_mask = torch.triu(torch.ones(time_lengh, time_lengh), 1) == 1

    # Calculate simulated FCD CDF
    fcd_hist = np.ones([10000, n_set])
    fc_mat = torch.zeros(fc_edgenum, sub_num, time_lengh)
    batch_num = int(n_set / sub_num)
    fc_resid = torch.zeros(fc_edgenum, resid_num, time_lengh)

    for b in range(batch_num):
        bold_temp = bold_d[:, b * sub_num:(b + 1) * sub_num, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, n_frame)
        for i in range(0, time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_maskm]
            fc_mat[:, :, i] = torch.transpose(
                cor_temp.view(sub_num, fc_edgenum), 0, 1)

        for j in range(0, sub_num):
            fcd_temp = torch_corr(torch.transpose(fc_mat[:, j, :], 0, 1))
            fcd_hist_temp = np.histogram(
                fcd_temp[fcd_mask].cpu().numpy(), bins=10000, range=(-1., 1.))
            fcd_hist[:, j + b * sub_num] = fcd_hist_temp[0]

    if resid_num != 0:
        bold_temp = bold_d[:, batch_num * sub_num:n_set, :]
        bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, n_frame)
        for i in range(time_lengh):
            bold_fc = torch_corr(bold_tempm[:, i:i + window_size])
            cor_temp = bold_fc[fc_mask_resid]
            fc_resid[:, :, i] = torch.transpose(
                cor_temp.view(resid_num, fc_edgenum), 0, 1)

        for j in range(resid_num):
            fcd_temp = torch_corr(torch.transpose(fc_resid[:, j, :], 0, 1))
            fcd_hist_temp = np.histogram(
                fcd_temp[fcd_mask].cpu().numpy(), bins=10000, range=(-1., 1.))
            fcd_hist[:, j + sub_num * batch_num] = fcd_hist_temp[0]

    fcd_histcum = np.cumsum(fcd_hist, 0)

    # Calculate KS cost
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
    print('Time taken to calculate FCD cost:', fcd_elapsed, 'seconds')
    return FCD_weight * ks_cost


def set_search_range(parser, FC):
    """

    This function defines the searching range of parameters (i.e. wEE, wEI,
    sigma) based on the input empirical FC.
    -Input:
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
        -FC: #ROI x #ROI empirical FC
    -Output:
        -search_range: (3x#ROI) x 2 2D matrix. searching range of each
        parameter (wEE, followed by wEI, followed by sigma). The first
        column sets the lower bound of the searching range, the second
        column sets the upper bound of the searching range.

    """

    n_node = FC.shape[0]
    dim = n_node * 3 + 1  # total number of parameters
    fc_mask = np.ones([n_node, n_node], dtype=bool)
    fc_mask = np.triu(fc_mask, 1)
    fc_mean = np.mean(FC[fc_mask])
    wEI_max = float(parser['training']['wEI_a']) * fc_mean + float(
        parser['training']['wEI_b'])
    wEI_min = float(parser['training']['wEI_c']) * fc_mean + float(
        parser['training']['wEI_d'])
    if wEI_min < 0:
        wEI_min = 0
    if wEI_min > 3.5:
        wEI_min = 3.5
    if wEI_max < 3.2:
        wEI_max = 3.2
    if wEI_max > 5:
        wEI_max = 5
    wEE_max = float(parser['training']['wEE_a']) * fc_mean + float(
        parser['training']['wEE_b'])
    wEE_min = float(parser['training']['wEE_c']) * fc_mean + float(
        parser['training']['wEE_d'])
    if wEE_min < 0:
        wEE_min = 0
    if wEE_min > 7:
        wEE_min = 7
    if wEE_max < 1:
        wEE_max = 1
    if wEE_max > 10:
        wEE_max = 10
    search_range = np.zeros((dim, 2))
    # search range for w_EE
    search_range[0:n_node, :] = [wEE_min, wEE_max]
    # search range for w_EI
    search_range[n_node:n_node * 2, :] = [wEI_min, wEI_max]
    # search range for G
    search_range[n_node * 2, :] = [0, 3]
    # search range for sigma
    search_range[n_node * 2 + 1:dim, :] = [0.0005, 0.01]
    return search_range


def csv_matrix_read(filename):
    """

    Convert .csv matrix to numpy array
    Input:
        -filename: the path of a file, string
    Output:
        -out_array: a numpy version of same array

    """

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


def torch_corr(A):
    """

    Matrix self-correlation function used with GPU

    """

    Amean = torch.mean(A, 1)
    Ax = A - torch.transpose(Amean.repeat(A.shape[1], 1), 0, 1)
    Astd = torch.mean(Ax**2, 1)
    Amm = torch.mm(Ax, torch.transpose(Ax, 0, 1)) / A.shape[1]
    Aout = torch.sqrt(torch.ger(Astd, Astd))
    Acor = Amm / Aout
    return Acor


def torch_corr2(A, B):
    """

    Compute correlation between 2 matrices (column wise)

    """

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


def CBIG_combined_cost_train(parameter, n_dup, SC, FC, FCD, rE_min, rE_max,
                             FC_weights, FCD_weight, parser):
    """

    A wrapper function used for training model parameters.
    Each execution of this function corresponds to 1 CMAES
    epoch.
    Input:
        -parameter: #feature x #set candidate parameter matrix.
            Each row corresponds to a feature (for example, wEE of the n-th
            ROI or global constant G)
        -n_dup: number of repetition simulated for each set of
            candidate parameters
        -SC: path to the training set SC (.csv file). #ROI x #ROI
        -FC: path to the training set FC (.csv file). #ROI x #ROI
        -FCD: path to the training set FCD (.mat file).
            The file should be a structure with a filed named 'FCD'.
            FCD_train should be a 10000x1 vector
        -rE_min, rE_max: lower and upper bounds of allowed
            excitatory firing rate thresholds, simulated BOLD with firing
            rate higher or lower than the thresholds will be set to NaN and
            excluded
        -FC_weights: 2 x 1 array. The first value is FC correlation loss
            weight, the second value is FC L1 loss weight
        -FCD_weight: 1x1 scalar
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
    -Output:
        -total_cost: total training cost associated with input model parameters
        -fc_corr_cost: training set FC correlation cost
        -fc_L1_cost: training set FC L1 cost
        -fcd_cost: training set FCD cost
        -bold: simulated BOLD signals
        -r_E_all: excitatory firing rate
    """

    simulation_period = float(parser['system']['simulation_period'])
    dt_training = float(parser['system']['dt_training'])
    low_mem = int(parser['training']['low_mem'])

    # convert to cuda
    parameter = torch.from_numpy(parameter).type(torch.DoubleTensor).cuda()

    # load training data
    emp_fcd = sio.loadmat(FCD)
    emp_fcd = np.array(emp_fcd['FCD'])

    sc_mat_raw = csv_matrix_read(SC)
    sc_mat = sc_mat_raw * 0.02 / sc_mat_raw.max()
    sc_mat = torch.from_numpy(sc_mat).type(torch.DoubleTensor).cuda()

    emp_fc = csv_matrix_read(FC)
    emp_fc = torch.from_numpy(emp_fc).type(torch.DoubleTensor).cuda()

    # run Euler integration
    bold, S_E_all, S_I_all, r_E_all, _ = CBIG_mfm_multi_simulation(
        parameter,
        sc_mat,
        simulation_period,
        n_dup,
        dt_training,
        low_mem,
        parser,
        rE_lower_bound=rE_min,
        rE_upper_bound=rE_max)

    # compute FC correlation and L1 losses
    fc_corr_cost, fc_L1_cost = CBIG_FC_multi_simulation(
        emp_fc, bold, n_dup, FC_weights)

    # compute FCD KS loss
    fcd_cost = CBIG_FCD_multi_simulation(emp_fcd, bold, n_dup, FCD_weight,
                                         parser)

    # compute total loss
    total_cost = fc_corr_cost + fc_L1_cost + fcd_cost

    return total_cost, fc_corr_cost, fc_L1_cost, fcd_cost, bold, r_E_all


def CBIG_combined_cost_validation(parameter, n_dup, SC, FC, FCD, rE_min,
                                  rE_max, FC_weights, FCD_weight, parser):
    """

     A wrapper function used for validating model parameters.
     epoch.
     Input:
         -parameter: #feature x #set candidate parameter matrix.
             Each row corresponds to a feature (for example, wEE of the n-th
             ROI or global constant G)
         -n_dup: number of repetition simulated for each set of
             candidate parameters
         -SC: path to the training set SC (.csv file). #ROI x #ROI
         -FC: path to the training set FC (.csv file). #ROI x #ROI
         -FCD: path to the training set FCD (.mat file). The file should be a
             structure with a filed named 'FCD'.
             FCD_train should be a 10000x1 vector
         -rE_min, rE_max: lower and upper bounds of allowed
             excitatory firing rate thresholds,  simulated BOLD with firing
             rate higher or lower than the thresholds will be set to NaN and
             excluded
         -FC_weights: 2 x 1 array. The first value is FC correlation
             loss weight, the second value is FC L1 loss weight
         -FCD_weight: 1x1 scalar
         -parser: a parsed configuration file specifying model parameters
             (see example.ini)
     -Output:
         -total_cost: total validation cost associated with input model
             parameters
         -fc_corr_cost: validation set FC correlation cost
         -fc_L1_cost: validation set FC L1 cost
         -fcd_cost: validation set FCD cost
         -bold: simulated BOLD signals
         -S_E_all: temporal average of excitatory synaptic gating variable
         -S_I_all: temporal average of inhibitory synaptic gating variable
         -r_E_all: excitatory firing rate
     """

    simulation_period = float(parser['system']['simulation_period'])
    dt_validation = float(parser['system']['dt_validation'])
    low_mem = int(parser['validation']['low_mem'])

    # convert to cuda
    parameter = torch.from_numpy(parameter).type(torch.DoubleTensor).cuda()

    # load validation data
    emp_fcd = sio.loadmat(FCD)
    emp_fcd = np.array(emp_fcd['FCD'])

    sc_mat_raw = csv_matrix_read(SC)
    sc_mat = sc_mat_raw * 0.02 / sc_mat_raw.max()
    sc_mat = torch.from_numpy(sc_mat).type(torch.DoubleTensor).cuda()

    emp_fc = csv_matrix_read(FC)
    emp_fc = torch.from_numpy(emp_fc).type(torch.DoubleTensor).cuda()

    # run Euler integration
    bold, S_E, S_I, r_E, J_I = CBIG_mfm_multi_simulation(
        parameter,
        sc_mat,
        simulation_period,
        n_dup,
        dt_validation,
        low_mem,
        parser,
        rE_lower_bound=rE_min,
        rE_upper_bound=rE_max)

    # compute FC correlation and L1 losses
    corr_cost, L1_cost = CBIG_FC_multi_simulation(emp_fc, bold, n_dup,
                                                  FC_weights)

    # compute FCD KS loss
    fcd_cost = CBIG_FCD_multi_simulation(emp_fcd, bold, n_dup, FCD_weight,
                                         parser)

    # compute total loss
    total_cost = corr_cost + L1_cost + fcd_cost

    return total_cost, corr_cost, L1_cost, fcd_cost, bold, S_E, S_I, r_E


def CBIG_combined_cost_test(parameter, n_dup, SC, FC, FCD, rE_min, rE_max,
                            FC_weights, FCD_weight, parser):
    """

     A wrapper function used for testing model parameters.
     epoch.
     Input:
         -parameter: #feature x #set candidate parameter matrix.
             Each row corresponds to a feature (for example, wEE of the n-th
             ROI or global constant G)
         -n_dup: number of repetition simulated for each set of
             candidate parameters
         -SC: path to the training set SC (.csv file). #ROI x #ROI
         -FC: path to the training set FC (.csv file). #ROI x #ROI
         -FCD: path to the training set FCD (.mat file).
             The file should be a structure with a filed named 'FCD'.
             FCD_train should be a 10000x1 vector
         -rE_min, rE_max: lower and upper bounds of allowed
             excitatory firing rate thresholds, simulated BOLD with firing
             rate higher or lower than the thresholds will be set to NaN and
             excluded
         -FC_weights: 2 x 1 array. The first value is FC correlation loss
             weight, the second value is FC L1 loss weight
         -FCD_weight: 1x1 scalar
         -parser: a parsed configuration file specifying model parameters
              (see example.ini)
     -Output:
         -total_cost: total training cost associated with input model
             parameters
         -fc_corr_cost: training set FC correlation cost
         -fc_L1_cost: training set FC L1 cost
         -fcd_cost: training set FCD cost
         -bold: simulated BOLD signals
         -r_E_all: excitatory firing rate
         -J_I: feedback inhibtion strength, i.e., wIE
     """

    simulation_period = float(parser['system']['simulation_period'])
    dt_test = float(parser['system']['dt_test'])

    # load test data
    parameter = np.tile(parameter, [1, n_dup])
    parameter = torch.from_numpy(parameter).type(torch.DoubleTensor).cuda()

    emp_fcd = sio.loadmat(FCD)
    emp_fcd = np.array(emp_fcd['FCD'])

    sc_mat_raw = csv_matrix_read(SC)
    sc_mat = sc_mat_raw * 0.02 / sc_mat_raw.max()
    sc_mat = torch.from_numpy(sc_mat).type(torch.DoubleTensor).cuda()

    emp_fc = csv_matrix_read(FC)
    emp_fc = torch.from_numpy(emp_fc).type(torch.DoubleTensor).cuda()

    # run Euler integration
    bold_d, S_E, S_I, r_E, J_I = CBIG_mfm_single_simulation(
        parameter,
        sc_mat,
        simulation_period,
        dt_test,
        parser,
        rE_lower_bound=rE_min,
        rE_upper_bound=rE_max)

    # compute FC correlation and L1 losses
    corr_cost, L1_cost = CBIG_FC_single_simulation(emp_fc, bold_d, n_dup,
                                                   FC_weights)

    # compute FCD KS loss
    fcd_cost = CBIG_FCD_single_simulation(emp_fcd, bold_d, n_dup, FCD_weight,
                                          parser)

    # compute total loss
    total_cost = corr_cost + L1_cost + fcd_cost

    return total_cost, corr_cost, L1_cost, fcd_cost, bold_d, S_E, S_I, r_E, J_I
