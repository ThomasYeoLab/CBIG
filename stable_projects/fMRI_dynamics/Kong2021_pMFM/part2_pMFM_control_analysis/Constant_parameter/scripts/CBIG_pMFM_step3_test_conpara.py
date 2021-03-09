# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import torch
import CBIG_pMFM_basic_functions as fc
import warnings


def CBIG_mfm_test_desikan_main(gpu_index=0):
    torch.cuda.set_device(gpu_index)
    torch.cuda.manual_seed(1)

    # Create output folder
    input_path = '../output/step2_validation_results/'
    output_path = '../output/step3_test_results/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    n_set = 100
    n_dup = 10

    n_node = 68
    vali_sel_num = 10
    vali_sel = np.zeros((3 * n_node + 1 + 8, vali_sel_num))

    for i in range(1, 11):
        load_file = 'random_seed_' + str(i) + '.csv'
        load_path = os.path.join(input_path, load_file)
        xmin = fc.csv_matrix_read(load_path)
        index_mat = np.zeros((2, xmin.shape[1]))
        index_mat[0, :] = i
        index_mat[1, :] = np.arange(xmin.shape[1])
        xmin = np.concatenate((index_mat, xmin), axis=0)
        vali_index = np.argmin(xmin[7, :])
        vali_sel[:, i - 1] = xmin[:, vali_index]

    vali_sort_index = np.argsort(vali_sel[7, :])
    vali_sel = vali_sel[:, vali_sort_index]

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
