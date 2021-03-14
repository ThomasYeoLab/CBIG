# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import torch
import scipy.io as sio
import CBIG_pMFM_basic_functions_main as fc
import warnings


def CBIG_pMFM_generate_simualted_fc_fcd(gpu_index=0):
    '''
    This function is to generate the simulated fc and fcd based on test set
    The simulated fc and fcd are used in the analysis shown in the paper

    Args:
        gpu_index:      index of gpu used for optimization
    Returns:
        None
    '''

    # Setting random seed and GPU
    torch.cuda.set_device(gpu_index)
    torch.cuda.manual_seed(1)

    # Create output folder
    test_file = '../output/step3_test_results/test_all.csv'
    output_path = '../output/step4_MFM_simulated_data'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    n_set = 2000
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

    # Calculating simualted BOLD signal using MFM
    bold_d = fc.CBIG_mfm_single_simulation(parameter, sc_mat, 14.4)

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

    # Calculating CDF for simualted FCD matrices
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

    fcd_dir = os.path.join(output_path, 'FCD')
    if not os.path.isdir(fcd_dir):
        os.makedirs(fcd_dir)
    tc_dir = os.path.join(output_path, 'TC')
    if not os.path.isdir(tc_dir):
        os.makedirs(tc_dir)

    count = 1
    for i in range(n_set):
        print('Generating simualted TC and FCD number: ' + str(count))
        fcd = fcd_numpy[:, :, i]
        bold = bold_numpy[:, i, :]
        if (fcd == fcd).all():
            FCD = {'FCD_mat': fcd}
            sio.savemat(
                os.path.join(fcd_dir, 'FCD_' + str(count) + '.mat'), FCD)
            BOLD = {'TC': bold}
            sio.savemat(
                os.path.join(tc_dir, 'TC_' + str(count) + '.mat'), BOLD)
            count += 1
        if count > 1000:
            break


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    CBIG_pMFM_generate_simualted_fc_fcd()
