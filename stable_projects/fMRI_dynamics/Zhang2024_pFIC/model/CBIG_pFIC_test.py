# /usr/bin/env python
"""
Written by Shaoshi Zhang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import shutil
import sys
import warnings
import configparser
import torch
import numpy as np
import scipy.io as sio
import CBIG_pFIC_misc as misc


def CBIG_mfm_test_main(parser):
    """

    This function runs forward simulation on the parameters with the lowest
    validation cost (referred as solutions). If multiple sets of solutions
    are considered (for example, run top 10 solutions with the lowest
    validation cost instead of just the top 1 solution), then a similarity
    constraint is applied. That is, the correlation between 2 solutions
    cannot be higher than a certain threshold (default is 0.98). By default,
    forward simulation is run 1000 times for each solution (each with a
    different noise instantiation) and the results are averaged every 10 times
    (so there are 100 costs in the end).
    -Input:
        -parser: a parsed configuration file specifying model parameters
        (see example.ini)
    -Output:
        -TC.mat: simulated BOLD time courses
        -S_E.mat: temporal average of simulated excitatory synaptic
            gating variable
        -S_I.mat: temporal average of simulated inhibitory synaptic
            gating variable
        -r_E.mat: temporal average of simulated excitatory firing rate
        -solution.csv: test costs of each repetition (each column represents
            an averaged FC corr cost, FC L1 cost and FCD cost across 10
            simulations, by default)
        -test_all.csv: a compiled file containing the training, validation
            and test costs for each solution. Each column corresponds to a
            solution. The first row corresponds to the process-id,
            the second row is the epoch number, the 3rd to 6th rows are
            training costs, 7th to 10th are validation costs, and 11th to
            14th are the test costs (averaged across all repetitions).

    """

    input_path = parser['validation']['output_path']
    output_path = parser['test']['output_path']
    FCD = parser['test']['FCD_test']
    SC = parser['test']['SC_test']
    FC = parser['test']['FC_test']
    rE_min = float(parser['system']['rE_min'])
    rE_max = float(parser['system']['rE_max'])
    FC_corr_weight = float(parser['training']['FC_corr_weight'])
    FC_L1_weight = float(parser['training']['FC_L1_weight'])
    FCD_weight = float(parser['training']['FCD_KS_weight'])
    gpu_number = int(parser['system']['GPU_index'])
    if gpu_number != -1:
        torch.cuda.set_device(gpu_number)

    if not os.path.exists(output_path):
        os.makedirs(output_path)
        os.makedirs(output_path + 'simulation/')
    else:
        shutil.rmtree(output_path)
        os.makedirs(output_path)
        os.makedirs(output_path + 'simulation/')

    torch.cuda.manual_seed(1)

    n_ave = int(parser['test']['n_ave'])
    n_dup = int(parser['test']['n_dup'])
    n_set = int(n_dup / n_ave)

    SC_data = misc.csv_matrix_read(SC)
    n_node = SC_data.shape[0]
    vali_raw_all = np.zeros((3 * n_node + 1 + 10, 1))

    num_process = len(eval(parser['system']['random_seed']))
    for i in range(1, num_process + 1):
        load_file = parser['validation']['output_file'] + '_' + str(i) + '.csv'
        load_path = os.path.join(input_path, load_file)
        if not os.path.exists(load_path):
            continue
        xmin = misc.csv_matrix_read(load_path)
        index_mat = np.zeros((2, xmin.shape[1]))
        index_mat[0, :] = i
        index_mat[1, :] = np.arange(xmin.shape[1])
        xmin = np.concatenate((index_mat, xmin), axis=0)

        vali_raw_all = np.concatenate((vali_raw_all, xmin), axis=1)

    vali_raw_all = vali_raw_all[:, 1:]
    vali_index = np.argsort(vali_raw_all[9, :])  # sort the total cost
    vali_sort_all = vali_raw_all[:, vali_index]

    n_solution = int(parser['test']['n_solution'])
    similarity_threshold = float(parser['test']['similarity_threshold'])
    i = 0
    vali_sel = np.zeros((vali_raw_all.shape[0], n_solution))
    p = 0
    p_set = np.zeros(n_solution)
    while i < n_solution and p < vali_raw_all.shape[1]:
        corr_t = np.zeros(n_solution, dtype=bool)
        corr_tr = np.zeros((n_solution, 3))
        for j in range(n_solution):
            wEE_corr = np.corrcoef(vali_sel[10:10 + n_node, j:j + 1].T,
                                   vali_sort_all[10:10 + n_node, p:p + 1].T)
            wEI_corr = \
                np.corrcoef(vali_sel[10 + n_node:10 + 2 * n_node, j:j + 1].T,
                            vali_sort_all[10 + n_node:10 + 2 * n_node,
                            p:p + 1].T)
            sigma_corr = \
                np.corrcoef(
                    vali_sel[11 + 2 * n_node:11 + 3 * n_node, j:j + 1].T,
                    vali_sort_all[11 + 2 * n_node:11 + 3 * n_node, p:p + 1].T)
            corr_tr[j, 0] = wEE_corr[0, 1]
            corr_tr[j, 1] = wEI_corr[0, 1]
            corr_tr[j, 2] = sigma_corr[0, 1]
        for k in range(n_solution):
            corr_t[k] = (np.abs(corr_tr[k, :] > similarity_threshold)).all()
        if not corr_t.any():
            vali_sel[:, i] = vali_sort_all[:, p]
            p_set[i] = p
            i += 1
        p += 1

    result_save = np.zeros((3 * n_node + 1 + 14, n_solution))
    result_save[0:10, :] = vali_sel[0:10, :]
    result_save[14:, :] = vali_sel[10:, :]

    for j in range(n_solution):
        print('[Solution ' + str(j + 1) + '] Test starts ...')
        test_cost = np.zeros((4, n_set))
        arx = np.tile(vali_sel[10:, j:j + 1], [1, n_set])
        total_cost, corr_cost, L1_cost, fcd_cost, bold, S_E, S_I, r_E, w_IE = \
            misc.CBIG_combined_cost_test(arx, n_ave, SC, FC, FCD, rE_min,
                                         rE_max, [FC_corr_weight,
                                                  FC_L1_weight], FCD_weight,
                                         parser)
        test_cost[0, :] = corr_cost
        test_cost[1, :] = L1_cost
        test_cost[2, :] = fcd_cost
        test_cost[3, :] = total_cost

        test_file = os.path.join(output_path,
                                 'solution_' + str(j + 1) + '.csv')
        with open(test_file, 'a'):
            np.savetxt(test_file, test_cost, delimiter=',')

        result_save[10, j] = np.nanmean(test_cost[0, :])
        result_save[11, j] = np.nanmean(test_cost[1, :])
        result_save[12, j] = np.nanmean(test_cost[2, :])
        result_save[13, j] = np.nanmean(test_cost[3, :])

        bold = bold.cpu().numpy()
        BOLD = {'TC': bold}
        S_E_save = S_E.cpu().numpy()
        S_E = {'S_E': S_E_save}
        S_I_save = S_I.cpu().numpy()
        S_I = {'S_I': S_I_save}
        r_E_save = r_E.cpu().numpy()
        r_E = {'r_E': r_E_save}

        sio.savemat(output_path + 'simulation/TC_' + str(j + 1) + '.mat', BOLD)
        sio.savemat(output_path + 'simulation/S_E_' + str(j + 1) + '.mat', S_E)
        sio.savemat(output_path + 'simulation/S_I_' + str(j + 1) + '.mat', S_I)
        sio.savemat(output_path + 'simulation/r_E_' + str(j + 1) + '.mat', r_E)
        print('*******************************************')

        test_file_all = os.path.join(output_path, 'test_all.csv')
        with open(test_file_all, 'a'):
            np.savetxt(test_file_all, result_save, delimiter=',')


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    config = sys.argv[1]
    parser = configparser.ConfigParser()
    parser.read(config)
    CBIG_mfm_test_main(parser)
