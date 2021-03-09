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


def CBIG_mfm_validation_desikan_main(gpu_index=0):
    '''
    This function is to validate the estimated parameters of mean field
    model.
    The objective function is the summation of FC correlation cost and
    FCD KS statistics cost.

    Args:
        gpu_index:      index of gpu used for optimization
        input_path:     input directory for the optimized model parameters
        output_path:    output directory for saving validation results
    Returns:
        None
    '''

    torch.cuda.set_device(gpu_index)

    # Create output folder
    input_path = '../output/step1_training_results/'
    output_path = '../output/step2_validation_results/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    highest_order = 1
    myelin_data = fc.csv_matrix_read('../input/myelin.csv')
    myelin_data = myelin_data[:, 0]
    n_node = myelin_data.shape[0]
    amatrix = np.zeros((n_node, highest_order + 1))
    for i in range(highest_order + 1):
        amatrix[:, i] = myelin_data**(i)
    template_mat = amatrix

    n_trial = 10
    vali_dup = 20

    for i in range(1, 11):
        random_seed_cuda = i + 100
        torch.cuda.manual_seed(random_seed_cuda)

        load_file = ['random_seed_', str(i), '.csv']
        load_path = [input_path] + load_file
        xmin = fc.csv_matrix_read(''.join(load_path))
        x_mass = xmin[0:-3, :]

        result_save = np.zeros((6 + 3 * n_node + 1, xmin.shape[1]))
        result_save[0:3, :] = xmin[-3:, :]

        para_w = template_mat @ x_mass[0:highest_order + 1, :]
        para_I = template_mat @ x_mass[highest_order + 1:2 *
                                       (highest_order + 1), :]
        para_sigma = template_mat @ x_mass[2 * (highest_order + 1) +
                                           1:x_mass.shape[0], :]
        arx_mass = np.concatenate(
            (para_w, para_I,
             x_mass[2 * (highest_order + 1):2 * (highest_order + 1) + 1, :],
             para_sigma), 0)

        result_save[6:, :] = arx_mass

        for k in range(n_trial):
            in_para = arx_mass[:, 50 * k:50 * (k + 1)]
            vali_total, vali_corr, vali_ks = \
                fc.CBIG_combined_cost_validation(
                    in_para, vali_dup)
            result_save[3, 50 * k:50 * (k + 1)] = vali_corr
            result_save[4, 50 * k:50 * (k + 1)] = vali_ks
            result_save[5, 50 * k:50 * (k + 1)] = vali_total

        save_path = [output_path] + load_file
        np.savetxt(''.join(save_path), result_save, delimiter=',')


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=UserWarning)
    CBIG_mfm_validation_desikan_main()
