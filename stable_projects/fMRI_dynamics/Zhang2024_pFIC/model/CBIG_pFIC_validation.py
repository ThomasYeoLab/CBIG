# /usr/bin/env python
"""
Written by Shaoshi Zhang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import math
import os
import shutil
import sys
from multiprocessing import Pool
from multiprocessing import set_start_method
import configparser
import torch
import numpy as np
import scipy.io as sio
import CBIG_pFIC_misc as misc


def CBIG_mfm_validation_main(parser, process_id):
    """

    All parameters generated at each CMAES epoch are passed through forward
    simulations to get validation costs.
    -Input:
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
        -process_id: a number used as an identifier for a process
    -Output:
        -TC.mat: simulated BOLD time courses
        -S_E.mat: temporal average of simulated excitatory
            synaptic gating variable
        -S_I.mat: temporal average of simulated inhibitory
            synaptic gating variable
        -r_E.mat: temporal average of  simulated excitatory firing rate
        [NOTE]: To save space, only the last 50 epoch's parameter's
            TC, S_E, S_I, r_E are saved as .mat
        -validation.csv: validation cost associated with each set of parameters
            Each column corresponds to 1 CMAES epoch. The first 4 rows are
            training FC corr loss, training FC L1 loss, training FCD KS loss
            and training total loss. The 5th to 8th rows are validation
            FC corr loss, validation FC L1 loss, validation FCD KS loss and
            validation total loss. The rest of the rows correspond to #ROI
            regional wEE, #ROI regional wEI, G, and #ROI regional sigma.
    """

    input_path = parser['training']['output_path']
    input_file = parser['training']['output_file'] + '_' + str(
        process_id) + '.csv'
    output_path = parser['validation']['output_path']
    output_file = parser['validation']['output_file'] + '_' + str(
        process_id) + '.csv'
    FCD = parser['validation']['FCD_validation']
    SC = parser['validation']['SC_validation']
    FC = parser['validation']['FC_validation']
    rE_min = float(parser['system']['rE_min'])
    rE_max = float(parser['system']['rE_max'])
    FC_corr_weight = float(parser['training']['FC_corr_weight'])
    FC_L1_weight = float(parser['training']['FC_L1_weight'])
    FCD_weight = float(parser['training']['FCD_KS_weight'])
    myelin = parser['training']['myelin_training']
    rsfc_gradient = parser['training']['RSFC_gradient_training']
    gpu_number = int(parser['system']['GPU_index'])
    if gpu_number != -1:
        torch.cuda.set_device(gpu_number)

    # create dir
    if not os.path.exists(output_path + 'simulation/' + str(process_id) + '/'):
        os.makedirs(output_path + 'simulation/' + str(process_id) + '/')
    else:
        shutil.rmtree(output_path + 'simulation/' + str(process_id) + '/')
        os.makedirs(output_path + 'simulation/' + str(process_id) + '/')

    myelin_data = misc.csv_matrix_read(myelin)
    num_myelin_component = myelin_data.shape[1]
    gradient_data = misc.csv_matrix_read(rsfc_gradient)
    num_gradient_component = gradient_data.shape[1]
    n_node = myelin_data.shape[0]
    N_p = num_myelin_component + num_gradient_component + 1
    concat_mat = np.vstack((np.ones(n_node), myelin_data.T, gradient_data.T)).T

    num_epoch = int(parser['training']['num_epoch'])
    n_trial = math.ceil(num_epoch / 50)
    vali_dup = 20

    random_seed_cuda = 1
    torch.cuda.manual_seed(random_seed_cuda)

    load_file = [input_file]
    load_path = [input_path] + load_file
    xmin = misc.csv_matrix_read(''.join(load_path))
    x_nonzero = xmin[:, xmin[0, :] != 0]
    x_mass = x_nonzero[0:-4, :]
    result_save = np.zeros((8 + 3 * n_node + 1, x_nonzero.shape[1]))
    result_save[0:4, :] = x_nonzero[-4:, :]
    wEE = concat_mat @ x_mass[0:N_p, :]
    wEI = concat_mat @ x_mass[N_p:2 * N_p, :]
    sigma = concat_mat @ x_mass[2 * N_p + 1:3 * N_p + 1, :]
    arx_mass = np.concatenate(
        (wEE, wEI, x_mass[2 * N_p:2 * N_p + 1, :], sigma), 0)
    result_save[8:, :] = arx_mass

    for k in range(n_trial):
        in_para = arx_mass[:, 50 * k:50 * (k + 1)]
        print("[Process ID: " + str(process_id) + "] Validation epoch " +
              str(k * 50) + " to epoch " + str((k + 1) * 50 - 1))
        vali_total, vali_corr, vali_L1, vali_ks, bold, S_E, S_I, r_E = \
            misc.CBIG_combined_cost_validation(in_para, vali_dup, SC, FC,
                                               FCD, rE_min, rE_max,
                                               [FC_corr_weight,
                                                FC_L1_weight], FCD_weight,
                                               parser)
        result_save[4, 50 * k:50 * (k + 1)] = vali_corr
        result_save[5, 50 * k:50 * (k + 1)] = vali_L1
        result_save[6, 50 * k:50 * (k + 1)] = vali_ks
        result_save[7, 50 * k:50 * (k + 1)] = vali_total

    bold_save = bold.cpu().numpy()
    BOLD = {'TC': bold_save}
    S_E_save = S_E.cpu().numpy()
    S_E = {'S_E': S_E_save}
    S_I_save = S_I.cpu().numpy()
    S_I = {'S_I': S_I_save}
    r_E_save = r_E.cpu().numpy()
    r_E = {'r_E': r_E_save}

    sio.savemat(output_path + 'simulation/' + str(process_id) + '/TC.mat',
                BOLD)
    sio.savemat(output_path + 'simulation/' + str(process_id) + '/S_E.mat',
                S_E)
    sio.savemat(output_path + 'simulation/' + str(process_id) + '/S_I.mat',
                S_I)
    sio.savemat(output_path + 'simulation/' + str(process_id) + '/r_E.mat',
                r_E)

    save_path = output_path + output_file
    with open(save_path, 'a'):
        np.savetxt(save_path, result_save, delimiter=',')


if __name__ == '__main__':
    config = sys.argv[1]
    parser = configparser.ConfigParser()
    parser.read(config)
    num_thread = int(parser['system']['num_thread'])

    input_args = ()
    for i in range(num_thread):
        input_args_thread = [parser, i + 1]
        input_args = input_args + (input_args_thread, )

    set_start_method("spawn")
    p = Pool(num_thread)
    p.starmap(CBIG_mfm_validation_main, input_args)
