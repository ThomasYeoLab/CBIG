"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch
import numpy as np
import scipy.io as spio
import os


def data_process_predictor(mode, n_epoch, group_mats_path, param_path, n_param=100, group_list=None):
    """The data processing function for DELSSOME FC+FCD cost predictor

    Args:
        mode (str): ['train', 'val', 'test']
        n_epoch (int): how many epochs of data you need to be put into the deep learning model
        group_mats_path (str): the path of the mat file saving myelin, gradient, SC, FC and FCD-CDF
        param_path (str): the path of the parameter saving directory
        n_param (int, optional): number of parameters in one epoch. Defaults to 100.
        group_list (list or 1D array, optional): the group_nbr list. Defaults to None.

    """

    torch.set_default_dtype(torch.float64)

    grouped_mats_path = os.path.join(group_mats_path, f'{mode}.mat')
    grouped_mats = spio.loadmat(grouped_mats_path)

    sc_sets = np.array(grouped_mats['sc_groups'])  # [n_set, n_roi, n_roi]
    sc_sets = torch.as_tensor(sc_sets)
    fc_sets = np.array(grouped_mats['fc_groups'])  # [n_set, n_roi, n_roi]
    fc_sets = torch.as_tensor(fc_sets)
    fcd_cdf_sets = np.array(grouped_mats['fcd_cdf_groups'])  # [n_set, bins]
    fcd_cdf_sets = torch.as_tensor(fcd_cdf_sets)

    if group_list is None:
        n_set = sc_sets.shape[0]
        group_list = np.arange(0, n_set)
    else:
        group_list = np.array(group_list)
        sc_sets = sc_sets[group_list]
        fc_sets = fc_sets[group_list]
        fcd_cdf_sets = fcd_cdf_sets[group_list]

    save_param_root = os.path.join(param_path, f'{mode}')
    param_save_dir_list = []

    for group_nbr in group_list:
        param_save_dir = os.path.join(save_param_root, f'group{group_nbr}')
        # Compatibility for pFIC (Deco2014) and pMFM (Deco2013) checkpoints: accept either file name.
        if not os.path.exists(os.path.join(param_save_dir, 'final_state_pFIC.pth')) and \
                not os.path.exists(os.path.join(param_save_dir, 'final_state.pth')):
            raise Exception(f"Final state does not exist in {param_save_dir}.")
        param_save_dir_list.append(param_save_dir)

    return data_processing_for_predictor(sc_sets, n_epoch, n_param, param_save_dir_list, fc_sets, fcd_cdf_sets)


def data_processing_for_predictor(sc_sets, n_epoch, n_param, param_save_dir_list, fc_sets, fcd_cdf_sets):
    """The data processing function for predict-simulate_empirical-loss model

    Args:
        sc_sets (array-like): [n_set, n_roi, n_roi]
        n_epoch (int): parameter epochs
        n_param (int): parameter number in each epoch
        param_save_dir_list (list of str): the list of parent directories,
            e.g. ['./sub0', './sub1', './sub2', ...]. The length must match n_set.
        fc_sets (array-like): [n_set, n_roi, n_roi].
        fcd_cdf_sets (array-like): [n_set, bins]. Will do normalization in this function

    Returns:
        tensors: parameter_sets, loss_sets, sc_sets, fc_sets, fcd_cdf_sets
    """

    device = 'cpu'

    n_set = sc_sets.shape[0]
    assert n_set == len(param_save_dir_list)
    n_roi = sc_sets.shape[1]
    parameter_dim = 3 * n_roi + 1

    sc_sets = torch.as_tensor(sc_sets)
    # sc_groups: [n_set, n_roi, n_roi]
    fc_sets = torch.as_tensor(fc_sets)
    fcd_cdf_sets = torch.as_tensor(fcd_cdf_sets)
    fcd_cdf_sets = fcd_cdf_sets / fcd_cdf_sets[:, -1].unsqueeze(1)  # Normalize for fcd_cdf

    parameter_sets = torch.ones(n_set, n_epoch, n_param, parameter_dim) * float('nan')
    loss_sets = torch.ones(n_set, n_epoch, n_param, 3) * float('nan')

    for set_nbr in range(n_set):
        for epoch in range(n_epoch):
            param_save_path = os.path.join(param_save_dir_list[set_nbr], f'param_save_epoch{epoch}.pth')
            if not os.path.exists(param_save_path):
                raise Exception(f"Parameter save path {param_save_path} does not exist.")
            d = torch.load(param_save_path, map_location=torch.device(device), weights_only=True)
            parameter = d['parameter']
            valid_param_list = d['valid_param_list']

            parameter_sets[set_nbr, epoch] = parameter.T
            loss_sets[set_nbr, epoch, valid_param_list] = d['FC_FCD_loss']  # [n_valid_param, 3]

    # Reshape parameter_sets to [n_set, epochs * n_param, parameter_dim]
    parameter_sets = parameter_sets.view(n_set, n_epoch * n_param, parameter_dim)
    # Reshape loss_sets to [n_set, epochs * n_param, 3]
    loss_sets = loss_sets.view(n_set, n_epoch * n_param, 3)

    return parameter_sets, loss_sets, sc_sets, fc_sets, fcd_cdf_sets
