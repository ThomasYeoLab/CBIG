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


def CBIG_mfm_extrapolation_main(parser):
    """

    If there is only fMRI partial coverage in the original data
    (like Alprazolam dataset), this function generates simulated time
    courses based on the extrapolated parameters covering the entire cortex.
    The extrapolation is done by multiplying the optimized linear
    coefficients with whole-cortex myelin and RSFC gradient.
    See this function for more details:
    [NOTE] Only the first solution will be considered
    -Input:
        -parser: a parsed configuration file specifying model parameters
        (see example.ini)
    -Output:
        -EI.mat: E/I ratio computed using the extrapolated time courses.

    """

    input_path = parser['test']['output_path'] + 'extrapolated_parameter.csv'
    output_path = parser['test']['output_path']
    SC = parser['extrapolation']['SC_full']
    rE_min = float(parser['system']['rE_min'])
    rE_max = float(parser['system']['rE_max'])
    gpu_number = int(parser['system']['GPU_index'])
    if gpu_number != -1:
        torch.cuda.set_device(gpu_number)
    torch.cuda.manual_seed(1)

    if not os.path.exists(output_path + 'extrapolation/'):
        os.makedirs(output_path + 'extrapolation/')
    else:
        shutil.rmtree(output_path + 'extrapolation/')
        os.makedirs(output_path + 'extrapolation/')
    n_set = 1000

    solution = misc.csv_matrix_read(input_path)
    solution = np.tile(solution, (1, n_set))
    solution = torch.from_numpy(solution).type(torch.DoubleTensor).cuda()

    SC = misc.csv_matrix_read(SC)
    SC = 0.02 * SC / SC.max()
    SC = torch.from_numpy(SC).type(torch.DoubleTensor).cuda()

    simulation_period = float(parser['system']['simulation_period'])
    dt_extrapolation = float(parser['system']['dt_extrapolation'])
    print('Extrapolation starts ...')
    _, S_E, S_I, r_E, _ = misc.CBIG_mfm_single_simulation(
        solution, SC, simulation_period, dt_extrapolation, parser)
    S_E = S_E.cpu().numpy()
    S_I = S_I.cpu().numpy()
    r_E = r_E.cpu().numpy()

    selection = np.ones(r_E.shape[1], dtype=np.int8)
    for i in range(r_E.shape[1]):
        if (rE_min > r_E[:, i]).any() or (r_E[:, i] > rE_max).any():
            selection[i] = False
    S_E = S_E[:, selection]
    S_I = S_I[:, selection]
    EI = np.mean(S_E, 1) / np.mean(S_I, 1)
    EI = {'EI': EI}
    sio.savemat(output_path + 'extrapolation/EI.mat', EI)


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    config = sys.argv[1]
    parser = configparser.ConfigParser()
    parser.read(config)
    CBIG_mfm_extrapolation_main(parser)
