# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import torch
import CBIG_pMFM_basic_functions_main as fc
import warnings


def CBIG_mfm_test_desikan_main(gpu_index=0):
    '''
    This function is to implement the testing processes of mean field
    model.
    The objective function is the summation of FC correlation cost and
    FCD KS statistics cost.

    Args:
        gpu_index:      index of gpu used for optimization
    Returns:
        None
    '''

    # Setting random seed and GPU
    torch.cuda.set_device(gpu_index)
    torch.cuda.manual_seed(1)

    # Create output folder
    input_path = '../output/step2_validation_results/'
    output_path = '../output/step3_test_results/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # Setting hyper-parameters
    n_set = 100
    n_dup = 10

    n_node = 68
    vali_raw_all = np.zeros((3 * n_node + 1 + 8, 1))

    for i in range(1, 11):
        load_file = 'random_seed_' + str(i) + '.csv'
        load_path = os.path.join(input_path, load_file)
        xmin = fc.csv_matrix_read(load_path)
        index_mat = np.zeros((2, xmin.shape[1]))
        index_mat[0, :] = i
        index_mat[1, :] = np.arange(xmin.shape[1])
        xmin = np.concatenate((index_mat, xmin), axis=0)

        vali_raw_all = np.concatenate((vali_raw_all, xmin), axis=1)

    vali_raw_all = vali_raw_all[:, 1:]
    vali_index = np.argsort(vali_raw_all[7, :])
    vali_sort_all = vali_raw_all[:, vali_index]

    vali_sel_num = 10
    i = 0
    vali_sel = np.zeros((vali_raw_all.shape[0], vali_sel_num))
    p = 0
    p_set = np.zeros(vali_sel_num)
    while i < vali_sel_num and p < vali_raw_all.shape[1]:
        corr_t = np.zeros(vali_sel_num, dtype=bool)
        corr_tr = np.zeros((vali_sel_num, 3))
        for j in range(vali_sel_num):
            w_corr = np.corrcoef(vali_sel[8:8 + n_node, j:j + 1].T,
                                 vali_sort_all[8:8 + n_node, p:p + 1].T)
            i_corr = np.corrcoef(
                vali_sel[8 + n_node:8 + 2 * n_node, j:j + 1].T,
                vali_sort_all[8 + n_node:8 + 2 * n_node, p:p + 1].T)
            s_corr = np.corrcoef(vali_sel[9 + 2 * n_node:, j:j + 1].T,
                                 vali_sort_all[9 + 2 * n_node:, p:p + 1].T)
            corr_tr[j, 0] = w_corr[0, 1]
            corr_tr[j, 1] = i_corr[0, 1]
            corr_tr[j, 2] = s_corr[0, 1]

        for k in range(vali_sel_num):
            corr_t[k] = (corr_tr[k, :] > 0.98).all()

        if not corr_t.any():
            vali_sel[:, i] = vali_sort_all[:, p]
            p_set[i] = p
            i += 1
        p += 1

    result_save = np.zeros((3 * n_node + 1 + 11, vali_sel_num))
    result_save[0:8, :] = vali_sel[0:8, :]
    result_save[11:, :] = vali_sel[8:, :]

    for j in range(vali_sel_num):
        test_cost = np.zeros((3, n_set * 10))
        for k in range(10):
            arx = np.tile(vali_sel[8:, j:j + 1], [1, n_set])
            total_cost, fc_cost, fcd_cost = fc.CBIG_combined_cost_test(
                arx, n_dup)
            test_cost[0, n_set * k:n_set * (k + 1)] = fc_cost
            test_cost[1, n_set * k:n_set * (k + 1)] = fcd_cost
            test_cost[2, n_set * k:n_set * (k + 1)] = total_cost
            test_file = os.path.join(output_path,
                                     'test_num_' + str(j + 1) + '.csv')
            np.savetxt(test_file, test_cost, delimiter=',')

        result_save[8, j] = np.nanmean(test_cost[0, :])
        result_save[9, j] = np.nanmean(test_cost[1, :])
        result_save[10, j] = np.nanmean(test_cost[2, :])
        print('****************  finish top ' + str(j + 1) +
              ' test  ****************')

        test_file_all = os.path.join(output_path, 'test_all.csv')
        np.savetxt(test_file_all, result_save, delimiter=',')


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=UserWarning)
    CBIG_mfm_test_desikan_main()
