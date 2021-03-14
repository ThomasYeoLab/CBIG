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


def CBIG_mfm_test_desikan_main_different_window(gpu_index=0,
                                                window_size='low'):
    '''
    This function is to compute the KS statistics for top 10 selected model
    parameters based on different sliding window length.

    Args:
        gpu_index:      index of gpu used for optimization
        window_size:    indicator of sliding window size:
                        low: 43
                        high: 125
    Returns:
        None
    '''
    input_path = '../../../part1_pMFM_main/output/' +\
                 'step3_test_results/test_all.csv'
    output_path = '../output/' + window_size
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    torch.cuda.set_device(gpu_index)
    torch.cuda.manual_seed(1)

    result_all = fc.csv_matrix_read(input_path)
    test_num = result_all.shape[1]

    n_set = 100
    n_dup = 10

    result_save = np.zeros((3, test_num))

    for j in range(0, 10):
        test_cost = np.zeros((3, 1000))
        for k in range(int(1000 / n_set)):
            arx = np.tile(result_all[11:, j:j + 1], [1, n_set])
            total_cost, fc_cost, fcd_cost = \
                fc.CBIG_combined_cost_test_differwin(arx, n_dup, window_size)
            test_cost[0, n_set * k:n_set * (k + 1)] = fc_cost
            test_cost[1, n_set * k:n_set * (k + 1)] = fcd_cost
            test_cost[2, n_set * k:n_set * (k + 1)] = total_cost
            test_file = os.path.join(output_path,
                                     'test_num_' + str(j + 1) + '.csv')
            np.savetxt(test_file, test_cost, delimiter=',')

        result_save[0, j] = np.nanmean(test_cost[0, :])
        result_save[1, j] = np.nanmean(test_cost[1, :])
        result_save[2, j] = np.nanmean(test_cost[2, :])
        print('**************  finish top ' + str(j + 1) +
              ' test  ***************')

        test_file_all = os.path.join(output_path, 'test_all.csv')
        np.savetxt(test_file_all, result_save, delimiter=',')


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    CBIG_mfm_test_desikan_main_different_window(window_size='low')
    CBIG_mfm_test_desikan_main_different_window(window_size='high')
