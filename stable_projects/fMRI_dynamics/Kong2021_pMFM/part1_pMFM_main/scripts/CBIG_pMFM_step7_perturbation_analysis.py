# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import torch
import time
import math
import scipy.io as sio
import CBIG_pMFM_basic_functions_main as fc
import warnings


def torch_max(A, B):
    if A.shape != B.shape:
        raise ValueError('Dimension mismatch.')
    Am = torch.unsqueeze(A, dim=len(A.shape))
    Bm = torch.unsqueeze(B, dim=len(B.shape))
    C = torch.cat((Am, Bm), dim=len(A.shape))
    o = torch.max(C, dim=len(A.shape))
    return o[0]


def torch_min(A, B):
    if A.shape != B.shape:
        raise ValueError('Dimension mismatch.')
    Am = torch.unsqueeze(A, dim=len(A.shape))
    Bm = torch.unsqueeze(B, dim=len(B.shape))
    C = torch.cat((Am, Bm), dim=len(A.shape))
    o = torch.min(C, dim=len(A.shape))
    return o[0]


def CBIG_mfm_original_simulation(parameter,
                                 sc_mat,
                                 t_epochlong,
                                 noise,
                                 d_t=0.01):
    '''
    Function used to generate the simulated BOLD signal using mean field
    model and hemodynamic model
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
    s_max = torch.zeros((n_nodes, n_set))
    s_min = torch.ones((n_nodes, n_set))
    cut_index = int(t_pre / t_bold)

    # Warm up
    start = time.time()
    for i in range(1000):
        d_y = fc.CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        y_t = y_t + d_y * d_t + w_coef * noise[:, :, i] * math.sqrt(d_t)

    # Main body: calculation
    for i in range(n_samples):
        d_y = fc.CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        random_num = noise[:, :, i + 1000]
        y_t = y_t + d_y * d_t + w_coef * random_num * math.sqrt(d_t)
        s_max = torch_max(y_t, s_max)
        s_min = torch_min(y_t, s_min)
        d_f = fc.CBIG_mfm_rfMRI_BW_ode(y_t, f_mat)
        f_mat = f_mat + d_f * d_t
        z_t, f_t, v_t, q_t = torch.chunk(f_mat, 4, dim=2)
        y_bold_temp = 100 / p_costant * v_0 * (
            k_1 * (1 - q_t) + k_2 * (1 - q_t / v_t) + k_3 * (1 - v_t))
        y_bold[:, :, count] = y_bold_temp[:, :, 0]
        count = count + ((i + 1) % (t_bold / d_t) == 0) * 1
    elapsed = time.time() - start
    print('The time used for calculating simulated BOLD signal is: ', elapsed)

    # Downsampling
    bold_d = y_bold[:, :, cut_index + 1:y_bold.shape[2]]

    return bold_d, s_max, s_min


def CBIG_mfm_perturbation_simulation(parameter,
                                     sc_mat,
                                     t_epochlong,
                                     noise,
                                     node_mask,
                                     index,
                                     svalue,
                                     d_t=0.01):
    '''
    Function used to generate the simulated BOLD signal using mean field
    model and hemodynamic model
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

    range_list = np.load('../output/perturbation_simulation/range_list.npy')
    range_index = range_list[index]
    start_point = range_index[0] + 60

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
    cut_index = int(t_pre / t_bold)

    # Warm up
    start = time.time()
    for i in range(1000):
        d_y = fc.CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        y_t = y_t + d_y * d_t + w_coef * noise[:, :, i] * math.sqrt(d_t)

    smax = torch.max(svalue[:, 0]) * node_mask
    smin = torch.min(svalue[:, 1]) * node_mask
    sign = 0
    # Main body: calculation
    for i in range(n_samples):
        d_y = fc.CBIG_mfm_rfMRI_ode(y_t, parameter, sc_mat)
        random_num = noise[:, :, i + 1000]
        y_t = y_t + d_y * d_t + w_coef * random_num * math.sqrt(d_t)
        if count >= start_point + cut_index and count < start_point + \
                cut_index + 1:
            y_t_masked = y_t * node_mask
            if torch.sum(y_t_masked) != torch.sum(y_t_masked):
                break
            if sign == 0 and torch.sum(abs(y_t_masked - smax)) <= torch.sum(
                    abs(y_t_masked - smin)):
                sign = -1

                def y_func(a):
                    return a
            elif sign == 0 and torch.sum(abs(y_t_masked - smax)) > torch.sum(
                    abs(y_t_masked - smin)):
                sign = 1

                def y_func(a):
                    return smax - a

            y_t = y_t + 0.8 * y_func(y_t_masked) * sign
        d_f = fc.CBIG_mfm_rfMRI_BW_ode(y_t, f_mat)
        f_mat = f_mat + d_f * d_t
        z_t, f_t, v_t, q_t = torch.chunk(f_mat, 4, dim=2)
        y_bold_temp = 100 / p_costant * v_0 * (
            k_1 * (1 - q_t) + k_2 * (1 - q_t / v_t) + k_3 * (1 - v_t))
        y_bold[:, :, count] = y_bold_temp[:, :, 0]
        count = count + ((i + 1) % (t_bold / d_t) == 0) * 1
    elapsed = time.time() - start
    print('The time used for calculating simulated BOLD signal is: ', elapsed)

    # Downsampling
    bold_d = y_bold[:, :, cut_index + 1:y_bold.shape[2]]

    return bold_d


def CBIG_pMFM_generate_simulated_original_data(gpu_index=0):
    torch.cuda.set_device(gpu_index)

    test_file = '../output/step3_test_results/test_all.csv'
    output_path = '../output/step7_perturbation_simulation/original'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    n_set = 200
    result_all = fc.csv_matrix_read(test_file)
    parameter = result_all[11:, 0]
    parameter = np.tile(parameter, [n_set, 1]).T
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    # Load data
    sc_mat_raw = fc.csv_matrix_read('../input/sc_test.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    count = 1

    for ti in range(10):
        print('Starting ' + str(ti))
        torch.cuda.manual_seed(ti)

        noise = torch.randn(68, n_set, 99402)

        # Calculating simulated BOLD signal using MFM
        bold_d, s_max, s_min = CBIG_mfm_original_simulation(
            parameter, sc_mat, 14.4, noise)

        # Initializing the FC and FCD masks
        n_set = bold_d.shape[1]
        n_nodes = bold_d.shape[0]
        window_size = 83
        time_length = 1200 - window_size + 1
        sub_num = 10
        fc_edgenum = int(n_nodes * (n_nodes - 1) / 2)
        fc_mask = torch.triu(torch.ones(n_nodes, n_nodes), 1) == 1
        fc_maskm = torch.zeros(n_nodes * sub_num,
                               n_nodes * sub_num).type(torch.cuda.ByteTensor)

        for i in range(sub_num):
            fc_maskm[n_nodes * i:n_nodes * (i + 1), n_nodes * i:n_nodes *
                     (i + 1)] = fc_mask

        # Calculating simulated FCD matrices
        fcd_all = torch.ones(time_length, time_length, n_set).cpu()
        fc_mat = torch.zeros(fc_edgenum, sub_num, time_length)
        batch_num = int(n_set / sub_num)

        for b in range(batch_num):
            bold_temp = bold_d[:, b * sub_num:(b + 1) * sub_num, :]
            bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, 1200)
            for i in range(0, time_length):
                bold_fc = fc.torch_corr(bold_tempm[:, i:i + window_size])
                cor_temp = bold_fc[fc_maskm]
                fc_mat[:, :, i] = torch.transpose(
                    cor_temp.view(sub_num, fc_edgenum), 0, 1)

            for j in range(0, sub_num):
                fcd_all[:, :, j + b * sub_num] = fc.torch_corr(
                    torch.transpose(fc_mat[:, j, :], 0, 1))

        bold_numpy = bold_d.cpu().numpy()
        fcd_numpy = fcd_all.numpy()
        noise_numpy = noise.cpu().numpy()
        smax_numpy = s_max.cpu().numpy()
        smin_numpy = s_min.cpu().numpy()

        # Save out simulated data
        fcd_dir = os.path.join(output_path, 'FCD')
        if not os.path.isdir(fcd_dir):
            os.makedirs(fcd_dir)
        tc_dir = os.path.join(output_path, 'TC')
        if not os.path.isdir(tc_dir):
            os.makedirs(tc_dir)
        noise_dir = os.path.join(output_path, 'Noise')
        if not os.path.isdir(noise_dir):
            os.makedirs(noise_dir)
        svalue_dir = os.path.join(output_path, 'Svalue')
        if not os.path.isdir(svalue_dir):
            os.makedirs(svalue_dir)

        for i in range(n_set):
            print('Generating simulated TC and FCD number: ' + str(count))
            fcd_save = fcd_numpy[:, :, i]
            bold_save = bold_numpy[:, i, :]
            noise_save = noise_numpy[:, i, :]
            svalue_save = np.zeros((n_nodes, 2))
            svalue_save[:, 0] = smax_numpy[:, i]
            svalue_save[:, 1] = smin_numpy[:, i]

            if (fcd_save == fcd_save).all():
                np.save(
                    os.path.join(fcd_dir, 'FCD_' + str(count) + '.npy'),
                    fcd_save)
                np.save(
                    os.path.join(tc_dir, 'TC_' + str(count) + '.npy'),
                    bold_save)
                np.save(
                    os.path.join(noise_dir, 'Noise_' + str(count) + '.npy'),
                    noise_save)
                np.save(
                    os.path.join(svalue_dir, 'Svalue_' + str(count) + '.npy'),
                    svalue_save)

                count += 1
            if count > 1000:
                break
        if count > 1000:
            break

        torch.cuda.empty_cache()


def CBIG_pMFM_determine_time_range():
    index_list = []
    range_list = []

    for index in range(1, 1001):

        FCD_mat = np.load(
            '../output/step7_perturbation_simulation/original/FCD/FCD_' +
            str(index) + '.npy')
        FCD_mean = np.mean(FCD_mat, 1)

        fcd_low = 1 * (FCD_mean < 0.6)
        len_count = 0
        max_count = 0
        range_index = np.array([0, 0])
        temp_start = 0
        for i in range(1, fcd_low.shape[0]):
            if fcd_low[i] == 1:
                len_count += 1
                if fcd_low[i - 1] == 0:
                    temp_start = i
            elif fcd_low[i] == 0 and fcd_low[i - 1] == 1:
                if max_count < len_count:
                    max_count = len_count
                    range_index[1] = i
                    range_index[0] = temp_start
                len_count = 0
        if max_count < len_count:
            max_count = len_count
            range_index[1] = i
            range_index[0] = temp_start

        # Only when the stable states last for more than 200 time step,
        # the FCD can be used in the perturbation experiment
        if max_count >= 200:
            index_list.append(index)
            range_list.append(range_index)

    # the index_list contains the indexes of FCD which can add perturbation
    # the range_list contains the perturbation injection time point
    np.save('../output/step7_perturbation_simulation/index_list.npy',
            index_list)
    np.save('../output/step7_perturbation_simulation/range_list.npy',
            range_list)


def CBIG_pMFM_generate_perturbed_FCD(gpu_index=0,
                                     region_num=5,
                                     region_indi='top'):
    test_file = '../output/step3_test_results/test_all.csv'
    output_path = '../output/step7_perturbation_simulation/' + \
                  region_indi + str(
                      region_num) + '_regions'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    torch.cuda.set_device(gpu_index)

    n_set = 1
    result_all = fc.csv_matrix_read(test_file)
    parameter = result_all[11:, 0]
    parameter = np.tile(parameter, [n_set, 1]).T
    parameter = torch.from_numpy(parameter).type(torch.FloatTensor).cuda()

    # Load data
    emp_fcd = sio.loadmat('../input/fcd_test.mat')
    emp_fcd = np.array(emp_fcd['test_aveM'])

    sc_mat_raw = fc.csv_matrix_read('../input/sc_test.csv')
    sc_mat = sc_mat_raw / sc_mat_raw.max() * 0.2
    sc_mat = torch.from_numpy(sc_mat).type(torch.FloatTensor).cuda()

    emp_fc = fc.csv_matrix_read('../input/fc_test.csv')
    emp_fc = torch.from_numpy(emp_fc).type(torch.FloatTensor).cuda()

    sim_grad_corr = sio.loadmat('../input/ROI_sim')
    sim_grad_corr = np.array(sim_grad_corr['SWSTD_FCD_sim'])
    sim_grad_corrM = np.tile(sim_grad_corr, [1, n_set])
    node_maskM = torch.from_numpy(sim_grad_corrM).cuda()
    sim_grad_corr_sort = np.sort(sim_grad_corr)
    if region_indi == 'top':
        node_mask = 1 * (node_maskM > sim_grad_corr_sort[-6]).type(
            torch.FloatTensor).cuda()
    else:
        node_mask = 1 * (node_maskM < sim_grad_corr_sort[5]).type(
            torch.FloatTensor).cuda()

    index_list = np.load(
        '../output/step7_perturbation_simulation/index_list.npy')

    fcd_dir = os.path.join(output_path, 'FCD')
    if not os.path.isdir(fcd_dir):
        os.makedirs(fcd_dir)
    tc_dir = os.path.join(output_path, 'TC')
    if not os.path.isdir(tc_dir):
        os.makedirs(tc_dir)

    for i in range(0, len(index_list)):
        index = index_list[i]
        print('Analyzing index ' + str(index))
        if os.path.isfile(os.path.join(fcd_dir, 'FCD_' + str(index) + '.npy')):
            continue

        noise_numpy = np.load(
            '../output/step7_perturbation_simulation/original/Noise'
            '/Noise_' + str(index) + '.npy')
        noise = torch.from_numpy(noise_numpy).cuda()
        noise = torch.unsqueeze(noise, dim=1)

        svalue_numpy = np.load(
            '../output/step7_perturbation_simulation/original/Svalue'
            '/Svalue_' + str(index) + '.npy')
        svalue = torch.from_numpy(svalue_numpy).type(torch.FloatTensor).cuda()

        # Calculating simulated BOLD signal using MFM
        bold_d = CBIG_mfm_perturbation_simulation(parameter, sc_mat, 14.4,
                                                  noise, node_mask, i, svalue)

        # Initializing the FC and FCD masks
        n_set = bold_d.shape[1]
        n_nodes = bold_d.shape[0]
        window_size = 83
        time_length = 1200 - window_size + 1
        sub_num = 1
        fc_edgenum = int(n_nodes * (n_nodes - 1) / 2)
        fc_mask = torch.triu(torch.ones(n_nodes, n_nodes), 1) == 1
        fc_maskm = torch.zeros(n_nodes * sub_num,
                               n_nodes * sub_num).type(torch.cuda.ByteTensor)

        for i in range(sub_num):
            fc_maskm[n_nodes * i:n_nodes * (i + 1), n_nodes * i:n_nodes *
                     (i + 1)] = fc_mask

        # Calculating CDF for simulated FCD matrices
        fcd_all = torch.ones(time_length, time_length, n_set).cpu()
        fc_mat = torch.zeros(fc_edgenum, sub_num, time_length)
        batch_num = int(n_set / sub_num)

        for b in range(batch_num):
            bold_temp = bold_d[:, b * sub_num:(b + 1) * sub_num, :]
            bold_tempm = bold_temp.transpose(0, 1).contiguous().view(-1, 1200)
            for i in range(0, time_length):
                bold_fc = fc.torch_corr(bold_tempm[:, i:i + window_size])
                cor_temp = bold_fc[fc_maskm]
                fc_mat[:, :, i] = torch.transpose(
                    cor_temp.view(sub_num, fc_edgenum), 0, 1)

            for j in range(0, sub_num):
                fcd_all[:, :, j + b * sub_num] = fc.torch_corr(
                    torch.transpose(fc_mat[:, j, :], 0, 1))

        bold_numpy = bold_d.cpu().numpy()
        fcd_numpy = fcd_all.numpy()

        fcd_save = fcd_numpy[:, :, 0]
        bold_save = bold_numpy[:, 0, :]

        np.save(os.path.join(fcd_dir, 'FCD_' + str(index) + '.npy'), fcd_save)
        np.save(os.path.join(tc_dir, 'TC_' + str(index) + '.npy'), bold_save)


def CBIG_pMFM_analysis_perturbed_FCD(region_num=5):
    index_list = np.load(
        '../output/step7_perturbation_simulation/index_list.npy')
    range_list = np.load(
        '../output/step7_perturbation_simulation/range_list.npy')

    origin_edges_all = np.array([])
    top_edges_all = np.array([])
    bottom_edges_all = np.array([])

    window_len = 200

    for i in range(0, index_list.shape[0]):
        index = index_list[i]

        fcd_origin = np.load(
            '../output/step7_perturbation_simulation/original/FCD/FCD_' +
            str(index) + '.npy')
        fcd_top = np.load('../output/step7_perturbation_simulation/top' +
                          str(region_num) + '/FCD/FCD_' + str(index) + '.npy')
        bold_top = np.load('../output/step7_perturbation_simulation/top' +
                           str(region_num) + '/TC/TC_' + str(index) + '.npy')
        fcd_bottom = np.load('../output/step7_perturbation_simulation/bottom' +
                             str(region_num) + '/FCD/FCD_' + str(index) +
                             '.npy')
        bold_bottom = np.load(
            '../output/step7_perturbation_simulation/bottom' +
            str(region_num) + '/TC/TC_' + str(index) + '.npy')

        if np.sum(bold_top[-1, :]) == 0 or np.isnan(np.sum(bold_top[-1, :])):
            continue
        if np.sum(bold_bottom[-1, :]) == 0 or np.isnan(
                np.sum(bold_bottom[-1, :])):
            continue

        range_index = range_list[i]
        perturb_start = range_index[0] + 18
        perturb_end = min(perturb_start + window_len, 1118)

        mat_origin = fcd_origin[perturb_start:perturb_end, perturb_start:
                                perturb_end]
        mat_top = fcd_top[perturb_start:perturb_end, perturb_start:perturb_end]
        mat_bottom = fcd_bottom[perturb_start:perturb_end, perturb_start:
                                perturb_end]

        origin_edges = np.mean(mat_origin, 1)
        top_edges = np.mean(mat_top, 1)
        bottom_edges = np.mean(mat_bottom, 1)

        origin_edges_all = np.concatenate((origin_edges_all,
                                           np.array([np.mean(origin_edges)])))
        top_edges_all = np.concatenate((top_edges_all,
                                        np.array([np.mean(top_edges)])))
        bottom_edges_all = np.concatenate((bottom_edges_all,
                                           np.array([np.mean(bottom_edges)])))

    output_dir = '../output/step7_perturbation_simulation/stats'
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    origin_data = {}
    origin_data['origin_edge'] = origin_edges_all

    top_data = {}
    top_data['top_edge'] = top_edges_all

    bottom_data = {}
    bottom_data['bottom_edge'] = bottom_edges_all

    sio.savemat(os.path.join(output_dir, 'origin_data.mat'), origin_data)
    sio.savemat(os.path.join(output_dir, 'top_data.mat'), top_data)
    sio.savemat(os.path.join(output_dir, 'bottom_data.mat'), bottom_data)


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=UserWarning)
    print('Start generating original siumulated data.')
    CBIG_pMFM_generate_simulated_original_data()
    print('Start determining perturbation starting time.')
    CBIG_pMFM_determine_time_range()
    print('Start generating perturbed simulated data')
    CBIG_pMFM_generate_perturbed_FCD(region_indi='top')
    CBIG_pMFM_generate_perturbed_FCD(region_indi='bottom')
    print('Start computing the final results')
    CBIG_pMFM_analysis_perturbed_FCD()
