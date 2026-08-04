"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import numpy as np
import torch
from models.CBIG_DELSSOME_pFIC import MfmModel


def get_EI_ratio(config, save_path, parameter, sc_euler, seed=None):
    """Compute E/I ratio from the input parameter and sc_euler

    Args:
        save_path (str): E/I ratio results saved path
        parameter (tensor): [parameter_dim, 1].
        param_dup (int): duplication times of parameter. The E/I ratio will be averaged across duplications
        sc_euler (tensor): structural connectivity matrix normalized to max = 0.02
        seed (int, optional): random seed for replication. Defaults to None.
    """
    # Set MFM hyper-parameters
    system_parameters = config['system']
    simulate_time = float(system_parameters['simulation_period'])
    burn_in_time = float(system_parameters['t_pre'])
    TR = float(system_parameters['TR'])
    param_dup = 3
    warm_up_t = int(system_parameters['warmup'])
    euler_dt = float(system_parameters['dt_test'])

    print(" -- Start simulating -- ")
    # Set random seed
    if seed is None:
        seed = np.random.randint(0, 1000000000000)
    torch.manual_seed(seed)

    # Apply parameter to simulate BOLD signals
    parameter_repeat = parameter.repeat(1, param_dup)  # [3*N+1, param_sets * param_dup]
    mfm_model = MfmModel(config, parameter_repeat, sc_euler, dt=euler_dt)
    bold_signal, valid_M_mask, s_e_ave, s_i_ave = mfm_model.CBIG_mfm_simulation(simulate_time=simulate_time,
                                                                                burn_in_time=burn_in_time,
                                                                                TR=TR,
                                                                                warm_up_t=warm_up_t,
                                                                                need_EI=True)
    print(f"Bold shape {bold_signal.shape}; S_E average shape {s_e_ave.shape}")
    # S_E: [n_roi, param_dup]

    # Calculate E/I ratio
    if valid_M_mask.any():
        s_e_this_param = s_e_ave[:, valid_M_mask]
        s_i_this_param = s_i_ave[:, valid_M_mask]
        ei_ratio = s_e_this_param / s_i_this_param  # [roi, valid_number]
        ei_ratio = torch.mean(ei_ratio, dim=1, keepdim=True)
        s_e_ave_ = torch.mean(s_e_this_param, dim=1, keepdim=True)  # [roi, 1]
        s_i_ave_ = torch.mean(s_i_this_param, dim=1, keepdim=True)
    else:
        print("No valid run.")
        return 1
    print("Start saving results...")
    save_dict = {
        'ei_ratio': ei_ratio,
        's_e_ave': s_e_ave_,
        's_i_ave': s_i_ave_,
        'parameter': parameter,
        'seed': seed,
        'time_series': bold_signal,
    }
    torch.save(save_dict, save_path)
    print("Successfully saved EI ratio.")
    print("-- E/I ratio Estimation Done. --")
    return 0
