"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import numpy as np
import torch
import torch.distributions as td
import os
import warnings

from models.CBIG_DELSSOME_model import MultiModalCostPredictor
from models.CBIG_DELSSOME_pFIC import MfmModel

from utils.CBIG_DELSSOME_pFIC_utils import parameterize_myelin_rsfc, CBIG_corr, FC_calculate, FCD_calculate, \
    all_loss_calculate_from_fc_fcd, read_predictor_constants


def train_help_function(config,
                        myelin,
                        RSFC_gradient,
                        sc_mat,
                        fc_emp,
                        fcd_cdf_emp,
                        save_param_dir,
                        num_epochs,
                        dl_epoch_range,
                        euler_epoch_range,
                        trials,
                        next_epoch,
                        seed=None):
    mfm_trainer = CMAES_Trainer(config=config,
                                myelin=myelin,
                                RSFC_gradient=RSFC_gradient,
                                sc_mat=sc_mat,
                                fc_emp=fc_emp,
                                fcd_cdf_emp=fcd_cdf_emp,
                                save_param_dir=save_param_dir,
                                num_epoch=num_epochs,
                                dl_epoch_range=dl_epoch_range,
                                euler_epoch_range=euler_epoch_range)
    if seed is None:
        trials = int(trials)
    else:
        trials = 1
    next_epoch = int(next_epoch)
    for i in range(trials):
        state = mfm_trainer.train_pFIC(seed=seed, next_epoch=next_epoch)
        # state: 0 success, 1 restart, 2 fail
        if state == 0:
            break
        elif state == 1:
            if i == trials - 1:
                print(f'Having tried for {trials} times, but still no available parameter.')
                break
        elif state == 2:
            print("CMA-ES broke during middle epochs. Terminate.")
            break
    return state


class CMAES_Trainer:

    def __init__(self, config, myelin, RSFC_gradient, sc_mat, fc_emp, fcd_cdf_emp, save_param_dir, num_epoch,
                 dl_epoch_range, euler_epoch_range):
        """Initialize CMA-ES trainer

        Args:
            config: config file content
            myelin (tensor): [ROIs, 1]
            RSFC_gradient (tensor): [ROIs, 1]
            sc_mat (tensor): structural connectivity. [ROIs, ROIs]
            fc_emp (tensor): empirical functional connectivity. [ROIs, ROIs]
            fcd_cdf_emp (tensor): empirical FCD CDF. Shape of [bins, 1], has been normalized.
            save_param_dir (str): the parameters saving directory
            num_epochs (int): total epochs
            dl_epoch_range (list): the epochs for DELSSOME CMA-ES
            euler_epoch_range (list): the epochs for Euler CMA-ES

        """

        self.device = 'cpu'

        self.config = config

        self.N = myelin.shape[0]
        self.myelin = myelin  # [N, 1]
        self.RSFC_gradient = RSFC_gradient  # [N, 1]
        self.concat_mat = torch.hstack((torch.ones_like(myelin), myelin, RSFC_gradient))  # [N, 3]
        self.pinv_concat_mat = torch.linalg.pinv(self.concat_mat)  # [3, N]

        self.sc_euler = torch.as_tensor(sc_mat / torch.max(sc_mat) * 0.02)
        # [N, N] for Euler integration
        self.sc_dl = self.sc_euler
        self.fc_euler = torch.as_tensor(fc_emp)
        self.fc_dl = self.fc_euler
        self.fcd_dl = torch.diff(fcd_cdf_emp.squeeze(), dim=0, prepend=torch.as_tensor([0]))  # [bins]
        self.fcd_cdf_euler = fcd_cdf_emp  # [bins, 1]

        self.parameter_dim = 3 * self.N + 1  # 3N + 1. The dimension of parameters in MFM model.
        self.d_param = 3  # The dimension of the parameter per region

        # System parameters
        system_parameters = config['system']
        self.simulate_time = float(system_parameters['simulation_period'])
        self.burn_in_time = float(system_parameters['t_pre'])
        self.TR = float(system_parameters['TR'])
        self.window_size = int(system_parameters['window_size'])
        assert self.N == int(system_parameters['n_ROI']), "Number of ROIs in config and myelin do not match."

        self.param_sets = 100
        # Every epoch the number of sampled parameter sets
        self.select_param_sets = 10
        # Number of selecting parameters to update CMA-ES from ${param_sets} parameter sets.
        self.min_select_param_sets = 1
        # If cannot get ${select_param_sets} parameters valid, the minimum number of valid parameters. Or CMA-ES break.
        self.param_dim = 10
        # The dimension of the parametrization coefficients
        self.param_dup = 3
        # duplicate for x times for each parameter in Euler integration to be more stable
        self.fcd_hist_bins = 10000
        # For reforming FCD matrix to probability distribution.
        self.warm_up_t = int(system_parameters['warmup'])
        # The duration of warm up loop in frames.
        self.dt = float(system_parameters['dt_training'])
        self.wEI_std_threshold = float(system_parameters['wEI_std_threshold'])
        if self.wEI_std_threshold == 0:
            self.wEI_std_threshold = None
        # If None, not consider, if float, the threshold will be applied to sampling procedure.

        # For wEI range constraint
        self.wEI_a = int(system_parameters['wEI_a'])
        self.wEI_b = int(system_parameters['wEI_b'])
        self.wEI_c = int(system_parameters['wEI_c'])
        self.wEI_d = int(system_parameters['wEI_d'])

        # For sigma range constraint
        self.sigma_scale = float(system_parameters['sigma_scale'])

        use_tqdm = int(system_parameters['use_tqdm'])
        if use_tqdm == 1:
            self.use_tqdm = True
        else:
            self.use_tqdm = False

        # Training and generating parameters
        self.epochs = int(num_epoch)
        self.dl_pfic_range = np.array(dl_epoch_range, dtype=int)  # 1D array
        self.euler_pfic_range = np.array(euler_epoch_range, dtype=int)
        self.pfic_range = np.concatenate((self.dl_pfic_range, self.euler_pfic_range))

        if len(self.dl_pfic_range) > 0:

            predictor_save_path, self.d_transformer = read_predictor_constants(config)
            self.dl_predictor = MultiModalCostPredictor(n_regions=self.N, d_param=self.d_param, d=self.d_transformer)
            self.dl_predictor.load_state_dict(
                torch.load(predictor_save_path, map_location=torch.device(self.device))['model_state_dict'])

        self.save_param_dir = save_param_dir
        if not os.path.exists(self.save_param_dir):
            os.makedirs(self.save_param_dir)
        print(f"""Successfully init CMA-ES trainer.
            The results will be saved under {self.save_param_dir}""")

    def get_parameters(self, param_10):
        """
        From the 10 parameterization parameters to get the 3N+1 parameters.
        The raw CMA-ES search vector (10-dim) is mapped to biologically constrained wEE/wEI/G/sigma
        via a linear projection using the pseudoinverse of the myelin/RSFC gradient concat matrix,
        following the MFM2013 parameterization.
        :param param_10: [10, param_sets]
        :return: parameters for 2014 Deco Model [3N+1, param_sets]
        """
        return parameterize_myelin_rsfc(self.myelin, self.RSFC_gradient, param_10)

    def get_wei_range(self):
        # Using mean FC to adjust wEI search range
        wei_max = self.wEI_a * torch.mean(self.fc_dl) + self.wEI_b
        wei_min = self.wEI_c * torch.mean(self.fc_dl) + self.wEI_d
        wei_min = min(max(wei_min, 0), 3.5)
        wei_max = min(max(wei_max, 3.2), 5)
        print(f"wEI parameter searching range: [{wei_min}, {wei_max}].")
        return wei_min, wei_max

    def mfm_model_loss_dl(self, parameter, epoch):
        """
        Using deep learning model to predict the total loss of MFM model
        :param parameter: [3*N+1, M]
        :param param_10: 10 parameterization parameters, just for saving.
        :param epoch: current epoch, just for naming the save file.
        :return: loss [select_param_sets(10),]; index [select_param_sets,]
        """
        self.dl_predictor.eval()
        with torch.no_grad():
            # Transpose: CMA-ES stores parameters as [features, batch]; predictor expects [batch, features]
            parameter = parameter.T  # [n_param, 3N+1]
            valid_param_list = torch.arange(0, parameter.shape[0], 1)
            count_valid = len(valid_param_list)
            select_param_sets = self.select_param_sets
            if count_valid >= self.select_param_sets:
                pass
            else:
                print(f'Valid parameter sets are not enough, only {count_valid} parameters!')
                if count_valid >= self.min_select_param_sets:
                    print('Choose min parameter sets instead.')
                    select_param_sets = count_valid
                else:
                    return None, None

            valid_parameter = parameter[valid_param_list]
            sc_this = self.sc_dl.unsqueeze(0).expand(valid_parameter.shape[0], -1, -1)
            fc_this = self.fc_dl.unsqueeze(0).expand(valid_parameter.shape[0], -1, -1)
            fcd_this = self.fcd_dl.unsqueeze(0).expand(valid_parameter.shape[0], -1)

            wEE = valid_parameter[:, 0:self.N]
            wEI = valid_parameter[:, self.N:2 * self.N]
            G = valid_parameter[:, 2 * self.N].unsqueeze(1).unsqueeze(1)  # [n_param, 1, 1]
            sigma = valid_parameter[:, 2 * self.N + 1:]

            # concatenate wEE, wEI, sigma as [n_param, n_roi, 3]
            parameter_sets = torch.stack((wEE, wEI, sigma), dim=-1)  # [n_param, n_roi, 3]
            sc_sets = G * sc_this  # [n_param, n_roi, n_roi]

            pred_loss = self.dl_predictor(parameter_sets, sc_sets, fc_this, fcd_this)  # [n_param, 3]
            total_loss = torch.sum(pred_loss, dim=1)  # [n_param]
            loss_sorted, index_sorted_in_valid = torch.sort(total_loss, descending=False)
            index_sorted = valid_param_list[index_sorted_in_valid]

            save_dict = {'parameter': parameter.T, 'valid_param_list': valid_param_list, 'FC_FCD_loss': pred_loss}
            torch.save(save_dict, os.path.join(self.save_param_dir, 'param_save_epoch' + str(epoch) + '.pth'))

            return loss_sorted[:select_param_sets], index_sorted[:select_param_sets]

    def mfm_model_loss_euler(self, parameter, epoch):
        """
        Calculate the total loss of MFM model
        :param parameter: [3*N+1, M]
        :param epoch: current epoch, just for naming the save file.
        :param need_regularization: 0 stands for pFIC, 1 stands for rFIC
        :return: loss [select_param_sets(10),]; index [select_param_sets,]
        """
        parameter_repeat = parameter.repeat(1, self.param_dup)  # [3*N+1, param_sets * param_dup]
        mfm_model = MfmModel(self.config, parameter_repeat, self.sc_euler, dt=self.dt)
        bold_signal, valid_M_mask = mfm_model.CBIG_mfm_simulation(simulate_time=self.simulate_time,
                                                                  burn_in_time=self.burn_in_time,
                                                                  TR=self.TR,
                                                                  warm_up_t=self.warm_up_t)
        # bold_signal: [N, M=param_sets * param_dup, t_for_bold]; valid_M_mask: [param_sets * param_dup]

        bold_signal = bold_signal.view(self.N, self.param_dup, self.param_sets, -1).transpose(1, 2)
        # [N, param_sets, param_dup, t_for_bold]
        valid_M_mask = valid_M_mask.view(self.param_dup, self.param_sets).T  # [param_sets, param_dup]

        valid_param_list = []  # record valid param index
        fc_sim = torch.zeros(self.param_sets, self.N, self.N)
        fcd_hist = torch.zeros(self.fcd_hist_bins, self.param_sets)
        count_valid = 0
        for i in range(self.param_sets):
            # for each set of parameter
            mask_this_param = valid_M_mask[i]  # [param_dup]
            if mask_this_param.any():
                valid_param_list.append(i)
                bold_this_param = bold_signal[:, i, mask_this_param, :]  # [N, 1/2/3/param_dup, t_for_bold]
                fc_this_param = FC_calculate(bold_this_param)
                fc_this_param = torch.mean(fc_this_param, dim=0)
                _, fcd_hist_this_param = FCD_calculate(bold_this_param, self.window_size)
                fcd_hist_this_param = torch.mean(fcd_hist_this_param, dim=1)

                fc_sim[count_valid] = fc_this_param
                fcd_hist[:, count_valid] = fcd_hist_this_param
                count_valid += 1

        select_param_sets = self.select_param_sets
        if count_valid >= self.select_param_sets:
            pass
        else:
            print(f'Valid parameter sets are not enough, only {count_valid} parameters!')
            if count_valid >= self.min_select_param_sets:
                print('Choose min parameter sets instead.')
                select_param_sets = count_valid
            else:
                return None, None

        fc_sim = fc_sim[:count_valid]
        fcd_hist = fcd_hist[:, :count_valid]
        total_loss, corr_loss, L1_loss, ks_loss = all_loss_calculate_from_fc_fcd(fc_sim, fcd_hist, self.fc_euler,
                                                                                 self.fcd_cdf_euler)  # [count_valid]
        FC_FCD_loss = torch.hstack((corr_loss.unsqueeze(1), L1_loss.unsqueeze(1), ks_loss.unsqueeze(1)))

        loss_sorted, index_sorted_in_valid = torch.sort(total_loss, descending=False)
        valid_param_list = torch.as_tensor(valid_param_list)
        index_sorted = valid_param_list[index_sorted_in_valid]

        save_dict = {'parameter': parameter, 'valid_param_list': valid_param_list, 'FC_FCD_loss': FC_FCD_loss}
        torch.save(save_dict, os.path.join(self.save_param_dir, 'param_save_epoch' + str(epoch) + '.pth'))

        return loss_sorted[:select_param_sets], index_sorted[:select_param_sets]

    def train_pFIC(self, seed=None, next_epoch=0):
        if next_epoch >= self.epochs:
            raise Exception('You do not need any more epoch.')

        N = self.N
        # Define search range. Maybe later it will depend on FC, so I do not put it in __init__
        search_range = torch.zeros(self.parameter_dim, 2)
        search_range[0:N, 0] = 1  # wee_min
        search_range[0:N, 1] = 10  # wee_max
        wei_min, wei_max = self.get_wei_range()
        search_range[N:2 * N, 0] = wei_min
        search_range[N:2 * N, 1] = wei_max
        # search range for G
        search_range[2 * N, 0] = 0
        search_range[2 * N, 1] = 3
        # search range for sigma
        search_range[2 * N + 1:, 0] = 0.0005 * self.sigma_scale
        search_range[2 * N + 1:, 1] = 0.01 * self.sigma_scale
        self.search_range = search_range

        if next_epoch == 0:
            if seed is None:
                seed = np.random.randint(0, 1000000000000)
            rand_ge = torch.manual_seed(seed)

            # initialization, k = 0
            m_k, sigma_k, cov_k, p_sigma_k, p_c_k = self._init_CMA_ES_pFIC(search_range)

        elif next_epoch in self.euler_pfic_range:
            previous_final_state_path = os.path.join(self.save_param_dir, 'final_state_pFIC.pth')
            if not os.path.exists(previous_final_state_path):
                raise Exception("Previous final state path doesn't exist.")
            final_dict = torch.load(previous_final_state_path, map_location=torch.device(self.device))
            seed = final_dict['seed']
            random_state = final_dict['random_state']
            m_k = final_dict['m']
            sigma_k = final_dict['sigma']
            cov_k = final_dict['cov']
            p_sigma_k = final_dict['p_sigma']
            p_c_k = final_dict['p_c']

            rand_ge = torch.manual_seed(seed)
            rand_ge.set_state(random_state)
            print("Successfully loaded previous parameters. Will start next epochs.")

        else:
            raise Exception("Argument next_epoch invalid .")

        for k in range(next_epoch, self.epochs):
            print("Epoch: [{}/{}]".format(k + 1, self.epochs))

            if k in self.pfic_range:
                epoch_res = self._train_one_epoch_pFIC(k, m_k, sigma_k, cov_k, p_sigma_k, p_c_k)
                if epoch_res is None:
                    # Return 1 (restart) if failure at the first epoch; return 2 (abort) for later epochs.
                    return 1 if k in self.dl_pfic_range or k == self.euler_pfic_range[0] else 2
                m_k, sigma_k, cov_k, p_sigma_k, p_c_k = epoch_res

        final_dict = {
            'seed': seed,
            'random_state': rand_ge.get_state(),
            'm': m_k,
            'sigma': sigma_k,
            'cov': cov_k,
            'p_sigma': p_sigma_k,
            'p_c': p_c_k,
            'epoch': k,
            'wEI_search_range': search_range[N:2 * N]
        }
        torch.save(final_dict, os.path.join(self.save_param_dir, 'final_state_pFIC.pth'))

        return 0

    def _train_one_epoch_pFIC(self, k, m_k, sigma_k, cov_k, p_sigma_k, p_c_k):
        param_10_k, parameter_k = self._sample_valid_parameters_pFIC(m_k, sigma_k**2 * cov_k, self.search_range)
        if param_10_k is None or parameter_k is None:
            print("Sampling failed!")
            return None

        if k in self.dl_pfic_range:
            loss_k, index_k = self.mfm_model_loss_dl(parameter_k, k)
        elif k in self.euler_pfic_range:
            loss_k, index_k = self.mfm_model_loss_euler(parameter_k, k)
        else:
            raise Exception("Invalid pFIC epoch. Break.")

        if loss_k is None:
            print("Abnormal loss encountered. End.")
            return None
        select_params = param_10_k[:, index_k]

        m_kp1, sigma_kp1, cov_kp1, p_sigma_kp1, p_c_kp1 = self._update_CMA_ES_pFIC(select_params, loss_k, m_k, sigma_k,
                                                                                   cov_k, p_sigma_k, p_c_k, k)
        return m_kp1, sigma_kp1, cov_kp1, p_sigma_kp1, p_c_kp1

    def _init_CMA_ES_pFIC(self, search_range):
        """
        Initialize from parameters' search range
        :return: m_0, cov_0, ... Their shape can be seen in self._update_CMA_ES
        """
        N = self.N

        # From 3N+1 parameters and concat_mat's pseudo inverse
        init_parameters = torch.rand(
            self.parameter_dim) * (search_range[:, 1] - search_range[:, 0]) + search_range[:, 0]  # [3*N+1]
        init_parameters = init_parameters.unsqueeze(1)  # [3*N+1, 1]
        start_point_wEE = torch.matmul(self.pinv_concat_mat, init_parameters[0:N]).squeeze()
        start_point_wEI = torch.matmul(self.pinv_concat_mat, init_parameters[N:2 * N]).squeeze()
        start_point_sigma = torch.matmul(self.pinv_concat_mat, init_parameters[2 * N + 1:]).squeeze()

        # Init m_0 for CMA-ES, just by experience
        m_0 = torch.zeros(self.param_dim)
        m_0[0:3] = start_point_wEE
        m_0[3:6] = start_point_wEI
        m_0[6] = init_parameters[2 * N]  # G
        m_0[7:] = start_point_sigma

        sigma_0 = 0.2
        p_sigma_0 = torch.zeros(self.param_dim, 1)
        p_c_0 = torch.zeros(self.param_dim, 1)
        V_ini = torch.eye(self.param_dim)
        Lambda_ini = torch.ones(self.param_dim)
        Lambda_ini[0:3] = start_point_wEE[0]
        Lambda_ini[3:6] = start_point_wEI[0]
        Lambda_ini[6] = 0.4
        Lambda_ini[7:] = 0.0005
        cov_0 = torch.matmul(V_ini, torch.matmul(torch.diag(Lambda_ini**2), V_ini.T))

        return m_0, sigma_0, cov_0, p_sigma_0, p_c_0

    def _update_CMA_ES_pFIC(self, select_params, loss_k, m_k, sigma_k, cov_k, p_sigma_k, p_c_k, k):
        """
        Update the parameters in CMA-ES algorithm (k -> k+1 (kp1)).
        Refer to wikipedia: https://en.wikipedia.org/wiki/CMA-ES
        :param select_params: parameters [param_dim, select_param_sets]
        :param loss_k: sorted loss by mfm_model_loss. [select_param_sets,]
        :param m_k: the means of multivariate gaussian, need updated. [param_dim,]
        :param sigma_k: the step sizes of CMA-ES, need updated. Double
        :param cov_k: the covariance matrix, need updated. [param_dim, param_dim]
        :param p_sigma_k: evolution path p_sigma
        :param p_c_k: evolution path p_c
        :param k: current iter
        :return: parameters in k+1 iter.
        """
        select_param_sets = select_params.shape[1]
        loss_inverse = 1 / loss_k
        weights = loss_inverse / torch.sum(loss_inverse)  # [select_param_sets, 1]
        m_kp1 = torch.matmul(select_params, weights)  # m_(k+1): [param_dim, 1]
        mueff = 1 / torch.sum(weights**2)  # mu_w

        # The evolution path p_sigma and p_c
        Lambda, V = torch.linalg.eigh(cov_k)  # eigen decomposition
        Lambda = torch.sqrt(Lambda)
        inv_sqrt_cov = torch.matmul(V, torch.matmul(torch.diag(Lambda**-1), V.T))  # C^(-1/2)

        c_sigma = (mueff + 2) / (self.param_dim + mueff + 5)
        c_c = (4 + mueff / self.param_dim) / (self.param_dim + 4 + 2 * mueff / self.param_dim)
        c_1 = 2 / ((self.param_dim + 1.3)**2 + mueff)
        c_mu = min(1 - c_1, 2 * (mueff - 2 + 1 / mueff) / ((self.param_dim + 2)**2 + mueff))
        d_sigma = 1 + 2 * max(0, torch.sqrt((mueff - 1) / (self.param_dim + 1)) - 1) + c_sigma
        expected_value = self.param_dim**0.5 * (1 - 1 / (4 * self.param_dim) + 1 / (21 * self.param_dim**2))

        p_sigma_kp1 = (1 - c_sigma) * p_sigma_k + torch.sqrt(c_sigma * (2 - c_sigma) * mueff) * torch.matmul(
            inv_sqrt_cov, (m_kp1 - m_k).unsqueeze(1) / sigma_k)
        indicator = (torch.linalg.norm(p_sigma_kp1) / torch.sqrt(1 - (1 - c_sigma)**(2 * k)) / expected_value
                     < (1.4 + 2 / (self.param_dim + 1))) * 1
        p_c_kp1 = (1 - c_c) * p_c_k + indicator * torch.sqrt(c_c *
                                                             (2 - c_c) * mueff) * (m_kp1 - m_k).unsqueeze(1) / sigma_k

        # Adapting covariance matrix C
        artmp = (1 / sigma_k) * (select_params - torch.tile(m_k, [select_param_sets, 1]).T)
        cov_kp1 = (1 - c_1 - c_mu) * cov_k + c_1 * (torch.matmul(p_c_kp1, p_c_kp1.T) + (1 - indicator) * c_c *
                                                    (2 - c_c) * cov_k) + c_mu * torch.matmul(
                                                        artmp, torch.matmul(torch.diag(weights), artmp.T))

        # Adapting step size
        sigma_kp1 = sigma_k * torch.exp((c_sigma / d_sigma) * (torch.linalg.norm(p_sigma_kp1) / expected_value - 1))
        return m_kp1, sigma_kp1, cov_kp1, p_sigma_kp1, p_c_kp1

    def _sample_valid_parameters_pFIC(self, mean, cov, search_range):
        # Reject sampled parameters that fail any of three biological constraints:
        # (a) parameter values outside the prescribed search range,
        # (b) wEI not negatively correlated with myelin or not positively correlated with RSFC gradient,
        # (c) wEI std below wEI_std_threshold (disabled when threshold is None), ensuring parameter diversity.
        multivariate_normal = td.MultivariateNormal(mean, cov)
        valid_count = 0
        total_count = 0
        total_threshold = 20000 * self.param_sets
        sampled_params = torch.zeros(self.param_dim, self.param_sets)
        sampled_parameters = torch.zeros(self.parameter_dim, self.param_sets)  # [3*N+1, param_sets]
        while valid_count < self.param_sets:
            sampled_params[:, valid_count] = multivariate_normal.sample()  # [10, param_sets]
            sampled_parameters[:, valid_count] = self.get_parameters(sampled_params[:, valid_count]).squeeze()
            wEI = sampled_parameters[self.N:2 * self.N, valid_count].unsqueeze(1)
            wEI_myelin_corr = torch.squeeze(CBIG_corr(wEI, self.myelin))
            wEI_rsfc_corr = torch.squeeze(CBIG_corr(wEI, self.RSFC_gradient))
            if self.wEI_std_threshold is None:
                if (sampled_parameters[:, valid_count] < search_range[:, 0]).any() \
                        or (sampled_parameters[:, valid_count] > search_range[:, 1]).any() \
                        or (wEI_myelin_corr > 0) or (wEI_rsfc_corr < 0):
                    valid_count -= 1
            else:
                if (sampled_parameters[:, valid_count] < search_range[:, 0]).any() \
                        or (sampled_parameters[:, valid_count] > search_range[:, 1]).any() \
                        or (wEI_myelin_corr > 0) or (wEI_rsfc_corr < 0) \
                        or (torch.std(wEI) < self.wEI_std_threshold):
                    valid_count -= 1

            valid_count += 1
            total_count += 1
            if total_count >= total_threshold:
                print(f"Not enough valid sampled parameters! Only sample {valid_count} parameters!")
                return None, None
        return sampled_params, sampled_parameters


class pFIC_Validator:

    def __init__(self, config, save_param_dir, val_save_dir, val_by_dl=False):
        """
        Args:
            config: config file content
            save_param_dir (str): the parameters saving directory from CMA-ES training
            val_save_dir (str): the saving directory for validation results
            val_by_dl (boolean): Defaults to False. whether to validate by DELSSOME model
        """

        self.device = 'cpu'

        self.config = config

        self.val_by_dl = val_by_dl

        # System parameters
        system_parameters = config['system']
        self.simulate_time = float(system_parameters['simulation_period'])
        self.burn_in_time = float(system_parameters['t_pre'])
        self.TR = float(system_parameters['TR'])
        self.window_size = int(system_parameters['window_size'])
        self.N = int(system_parameters['n_ROI'])
        self.parameters_dim = 3 * self.N + 1  # The dimension of parameters in MFM model.
        self.d_param = 3  # The dimension of the parameter per region
        self.param_sets = 100
        self.param_dup = 3
        self.fcd_hist_bins = 10000
        self.warm_up_t = int(system_parameters['warmup'])
        self.dt = float(system_parameters['dt_validation'])

        if self.val_by_dl:

            predictor_save_path, d_transformer = read_predictor_constants(config)
            self.dl_predictor = MultiModalCostPredictor(n_regions=self.N, d_param=self.d_param, d=d_transformer)
            self.dl_predictor.load_state_dict(
                torch.load(predictor_save_path, map_location=torch.device(self.device))['model_state_dict'])

        self.save_param_dir = save_param_dir
        self.val_save_dir = val_save_dir

        if not os.path.exists(self.val_save_dir):
            os.makedirs(self.val_save_dir)
        print(f"""pFIC validator has been successfully initialized.
            Results will be saved in {self.val_save_dir}""")

    def val_best_parameters_all_epochs(self, sc_euler, fc_emp, fcd_cdf_emp, num_epochs, seed=None):
        # Set random seed
        if seed is None:
            seed = np.random.randint(0, 1000000000000)
        torch.manual_seed(seed)

        parameter_all_epochs = torch.zeros(self.parameters_dim, num_epochs)
        for epoch in range(num_epochs):
            parameter_path = os.path.join(self.save_param_dir, f'param_save_epoch{epoch}.pth')
            if not os.path.exists(parameter_path):
                warnings.warn("Path doesn't exist.")
                return 1
            d = torch.load(parameter_path, map_location=torch.device(self.device))

            valid_param_list_train = d['valid_param_list']
            parameter = d['parameter']
            parameter = parameter[:, valid_param_list_train]
            if 'FC_FCD_loss' in d:
                total_loss = torch.sum(d['FC_FCD_loss'], dim=1)  # [n_param]
            else:
                raise Exception("Check the dictionary keys.")
            best_param_ind = torch.argmin(total_loss)
            parameter_all_epochs[:, epoch] = parameter[:, best_param_ind]

        if self.val_by_dl:
            losses = self._val_best_parameters_all_epochs_by_dl(parameter_all_epochs, sc_euler, fc_emp, fcd_cdf_emp,
                                                                num_epochs)
        else:
            losses = self._val_best_parameters_all_epochs_by_euler(parameter_all_epochs, sc_euler, fc_emp, fcd_cdf_emp,
                                                                   num_epochs)
        if losses is None:
            return 1
        else:
            _, corr_loss, L1_loss, ks_loss, valid_param_list = losses

        save_dict = {
            'parameter': parameter_all_epochs,
            'valid_param_list': valid_param_list,
            'corr_loss': corr_loss,
            'l1_loss': L1_loss,
            'ks_loss': ks_loss,
            'seed': seed
        }
        torch.save(save_dict, os.path.join(self.val_save_dir, 'best_params_all_epochs.pth'))
        return 0

    def _val_best_parameters_all_epochs_by_euler(self, parameter_all_epochs, sc_euler, fc_emp, fcd_cdf_emp, num_epochs):
        parameter_repeat = parameter_all_epochs.repeat(1, self.param_dup)  # [3*N+1, param_sets * param_dup]
        mfm_model = MfmModel(self.config, parameter_repeat, sc_euler, dt=self.dt)
        bold_signal, valid_M_mask = mfm_model.CBIG_mfm_simulation(simulate_time=self.simulate_time,
                                                                  burn_in_time=self.burn_in_time,
                                                                  TR=self.TR,
                                                                  warm_up_t=self.warm_up_t)

        bold_signal = bold_signal.view(self.N, self.param_dup, num_epochs, -1).transpose(1, 2)
        # [N, param_sets, param_dup, t_for_bold]
        valid_M_mask = valid_M_mask.view(self.param_dup, num_epochs).T  # [param_sets, param_dup]

        valid_param_list = []  # record valid param index
        fc_sim = torch.zeros(num_epochs, self.N, self.N)
        fcd_hist = torch.zeros(self.fcd_hist_bins, num_epochs)
        count_valid = 0
        for i in range(num_epochs):
            # for each set of parameter
            mask_this_param = valid_M_mask[i]  # [param_dup]
            if mask_this_param.any():
                valid_param_list.append(i)
                bold_this_param = bold_signal[:, i, mask_this_param, :]  # [N, 1/2/3/param_dup, t_for_bold]
                fc_this_param = FC_calculate(bold_this_param)
                fc_this_param = torch.mean(fc_this_param, dim=0)
                _, fcd_hist_this_param = FCD_calculate(bold_this_param, self.window_size)
                fcd_hist_this_param = torch.mean(fcd_hist_this_param, dim=1)

                fc_sim[count_valid] = fc_this_param
                fcd_hist[:, count_valid] = fcd_hist_this_param
                count_valid += 1
        fc_sim = fc_sim[:count_valid]
        fcd_hist = fcd_hist[:, :count_valid]
        total_loss, corr_loss, L1_loss, ks_loss = all_loss_calculate_from_fc_fcd(fc_sim, fcd_hist, fc_emp,
                                                                                 fcd_cdf_emp)  # [count_valid]
        valid_param_list = torch.as_tensor(valid_param_list)
        return total_loss, corr_loss, L1_loss, ks_loss, valid_param_list

    def _val_best_parameters_all_epochs_by_dl(self, parameter_all_epochs, sc_euler, fc_emp, fcd_cdf_emp, num_epochs):
        self.dl_predictor.eval()
        with torch.no_grad():
            parameter_all_epochs = parameter_all_epochs.T  # [param_sets, 3N+1]
            sc_this = sc_euler.unsqueeze(0).expand(parameter_all_epochs.shape[0], -1, -1)
            fc_this = fc_emp.unsqueeze(0).expand(parameter_all_epochs.shape[0], -1, -1)
            fcd_this = torch.diff(fcd_cdf_emp.squeeze(), dim=0, prepend=torch.as_tensor([0])).unsqueeze(0)
            # [1, 10000]
            fcd_this = fcd_this.expand(parameter_all_epochs.shape[0], -1)

            wEE = parameter_all_epochs[:, 0:self.N]
            wEI = parameter_all_epochs[:, self.N:2 * self.N]
            G = parameter_all_epochs[:, 2 * self.N].unsqueeze(1).unsqueeze(1)  # [n_param, 1, 1]
            sigma = parameter_all_epochs[:, 2 * self.N + 1:]

            # concatenate wEE, wEI, sigma as [n_param, n_roi, 3]
            parameter_sets = torch.stack((wEE, wEI, sigma), dim=-1)  # [n_param, n_roi, 3]
            sc_sets = G * sc_this  # [n_param, n_roi, n_roi]

            pred_loss = self.dl_predictor(parameter_sets, sc_sets, fc_this, fcd_this)  # [n_param, 3]
            total_loss = torch.sum(pred_loss, dim=1)  # [n_param]
            return total_loss, pred_loss[:, 0], pred_loss[:, 1], pred_loss[:, 2], torch.arange(0, num_epochs, 1)


class pFIC_Tester:

    def __init__(self, config, val_dirs, test_dir, trained_epochs):
        """
        Args:
            config: config file content
            val_dirs: directories of validation/train results. Mainly from multiple seeds.
            [NOTE]: when val_dirs are train directories, use function select_best_from_train
            when val_dirs are validation directories, use function select_best_from_val
            test_dir: diretory saving test results.
            trained_epochs (int): total epochs
        """

        self.device = 'cpu'

        self.config = config

        system_parameters = config['system']
        self.simulate_time = float(system_parameters['simulation_period'])
        self.burn_in_time = float(system_parameters['t_pre'])
        self.TR = float(system_parameters['TR'])
        self.window_size = int(system_parameters['window_size'])
        self.N = int(system_parameters['n_ROI'])
        self.parameters_dim = 3 * self.N + 1  # The dimension of parameters in MFM model.
        self.param_sets = 1  # select number of best parameters from val_dirs
        # [NOTE] here param_sets is not identical to the definition in Validator and Trainer
        self.param_dup = 3
        self.fcd_hist_bins = 10000
        self.warm_up_t = int(system_parameters['warmup'])
        self.dt = float(system_parameters['dt_test'])

        self.trained_epochs = trained_epochs  # The number of param files in one validation directory

        self.val_dirs = val_dirs
        self.test_dir = test_dir
        self.val_dirs_len = len(self.val_dirs)

        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)
        print(f"""pFIC testing process successfully initialized.
        Results will be saved under {self.test_dir}""")

    def test(self, sc_euler, fc_emp, fcd_cdf_emp, seed=None):
        """Apply best performed parameters in validation set to test set

        Args:
            sc_euler (tensor): structural connectivity normalized to max = 0.02. [ROIs, ROIs]
            fc_emp (tensor): empirical functional connectivity. [ROIs, ROIs]
            fcd_cdf_emp (tensor): empirical FCD CDF. [bins, 1]
            seed (int, optional): random seed for replication. Defaults to None.

        Returns:
            int: state.
        """

        # Set random seed
        if seed is None:
            seed = np.random.randint(0, 1000000000000)
        torch.manual_seed(seed)

        # Load validation results
        parameter_sets = torch.zeros(self.parameters_dim, self.val_dirs_len * self.trained_epochs)
        val_loss_sets = torch.ones(self.val_dirs_len * self.trained_epochs) * 3
        valid_val_dir_count = 0
        valid_val_epoch_count = 0
        for val_dir_i in range(self.val_dirs_len):
            val_dir = self.val_dirs[val_dir_i]
            if not os.path.exists(val_dir):
                print(f"{val_dir} cannot be found.")
                continue
            valid_val_dir_count += 1
            param_all_epochs_val_path = os.path.join(val_dir, 'best_params_all_epochs.pth')
            if os.path.exists(param_all_epochs_val_path):
                d = torch.load(param_all_epochs_val_path, map_location=torch.device(self.device))
                if len(d['valid_param_list']) == 0:
                    print("valid_param_list is empty.")
                    continue
                parameter_sets[:,
                               val_dir_i * self.trained_epochs:(val_dir_i + 1) * self.trained_epochs] = d['parameter']
                val_loss_sets[val_dir_i * self.trained_epochs:(val_dir_i + 1) *
                              self.trained_epochs][d['valid_param_list']] = d['corr_loss'] + d['l1_loss'] + d['ks_loss']
                valid_val_epoch_count += len(d['valid_param_list'])
            else:
                print(f"{param_all_epochs_val_path} cannot be found. Looking for best_param per epoch.")
                for epoch in range(self.trained_epochs):
                    param_val_path = os.path.join(val_dir, f'best_param{epoch}.pth')
                    if not os.path.exists(param_val_path):
                        print(f"{param_val_path} cannot be found.")
                        continue
                    valid_val_epoch_count += 1
                    d = torch.load(param_val_path, map_location=torch.device(self.device))
                    parameter_sets[:, val_dir_i * self.trained_epochs + epoch] = torch.squeeze(d['parameter'])
                    val_loss_sets[val_dir_i * self.trained_epochs +
                                  epoch] = d['corr_loss'] + d['l1_loss'] + d['ks_loss']
        if valid_val_dir_count == 0:
            print("No valid validated directories.")
            return 1
        if valid_val_epoch_count == 0:
            print("No valid epoch.")
            return 1

        # Record and sort all parameters and loss
        val_losses, sorted_index = torch.sort(val_loss_sets, descending=False)
        val_losses = val_losses[:self.param_sets]
        sorted_index = sorted_index[:self.param_sets]
        parameter = parameter_sets[:, sorted_index]

        # Apply the best performed ${param_sets} parameters to test set
        parameter_repeat = parameter.repeat(1, self.param_dup)  # [3*N+1, param_sets * param_dup]
        mfm_model = MfmModel(self.config, parameter_repeat, sc_euler, dt=self.dt)
        bold_signal, valid_M_mask = mfm_model.CBIG_mfm_simulation(simulate_time=self.simulate_time,
                                                                  burn_in_time=self.burn_in_time,
                                                                  TR=self.TR,
                                                                  warm_up_t=self.warm_up_t)

        bold_signal = bold_signal.view(self.N, self.param_dup, self.param_sets, -1).transpose(1, 2)
        # [N, param_sets, param_dup, t_for_bold]
        valid_M_mask = valid_M_mask.view(self.param_dup, self.param_sets).T  # [param_sets, param_dup]

        # Compute FC, FCD and then loss
        valid_param_list = []  # record valid param index
        fc_sim = torch.zeros(self.param_sets, self.N, self.N)
        fcd_hist = torch.zeros(self.fcd_hist_bins, self.param_sets)
        count_valid = 0
        for i in range(self.param_sets):
            # for each set of parameter
            mask_this_param = valid_M_mask[i]  # [param_dup]
            if mask_this_param.any():
                valid_param_list.append(i)
                bold_this_param = bold_signal[:, i, mask_this_param, :]  # [N, 1/2/3/param_dup, t_for_bold]
                fc_this_param = FC_calculate(bold_this_param)
                fc_this_param = torch.mean(fc_this_param, dim=0)
                _, fcd_hist_this_param = FCD_calculate(bold_this_param, self.window_size)
                fcd_hist_this_param = torch.mean(fcd_hist_this_param, dim=1)

                fc_sim[count_valid] = fc_this_param
                fcd_hist[:, count_valid] = fcd_hist_this_param
                count_valid += 1
        fc_sim = fc_sim[:count_valid]
        fcd_hist = fcd_hist[:, :count_valid]
        total_loss, corr_loss, L1_loss, ks_loss = all_loss_calculate_from_fc_fcd(fc_sim, fcd_hist, fc_emp,
                                                                                 fcd_cdf_emp)  # [count_valid]
        valid_param_list = torch.as_tensor(valid_param_list)

        save_dict = {
            'parameter': parameter,
            'val_losses': val_losses,
            'sorted_index': sorted_index,
            'valid_param_list': valid_param_list,
            'time_series': bold_signal,
            'corr_loss': corr_loss,
            'l1_loss': L1_loss,
            'ks_loss': ks_loss,
            'seed': seed
        }
        torch.save(save_dict, os.path.join(self.test_dir, 'test_results.pth'))
        return 0

    def select_best_from_val(self):
        """Select the best performed parameter in validation set, normally a step before E/I ratio calculation

        Returns:
            int: state.
        """

        parameter_sets = torch.zeros(self.parameters_dim, self.val_dirs_len * self.trained_epochs)
        val_loss_sets = torch.ones(self.val_dirs_len * self.trained_epochs, 4) * 3
        # [total, corr, L1, ks]
        valid_val_dir_count = 0
        valid_val_epoch_count = 0
        for val_dir_i in range(self.val_dirs_len):
            val_dir = self.val_dirs[val_dir_i]
            if not os.path.exists(val_dir):
                print(f"{val_dir} cannot be found.")
                continue
            valid_val_dir_count += 1
            param_all_epochs_val_path = os.path.join(val_dir, 'best_params_all_epochs.pth')
            if os.path.exists(param_all_epochs_val_path):
                d = torch.load(param_all_epochs_val_path, map_location=torch.device(self.device))
                if len(d['valid_param_list']) == 0:
                    print("valid_param_list is empty.")
                    continue
                parameter_sets[:,
                               val_dir_i * self.trained_epochs:(val_dir_i + 1) * self.trained_epochs] = d['parameter']
                val_loss_sets[val_dir_i * self.trained_epochs:(val_dir_i + 1) *
                              self.trained_epochs][d['valid_param_list'], 1] = d['corr_loss']
                val_loss_sets[val_dir_i * self.trained_epochs:(val_dir_i + 1) *
                              self.trained_epochs][d['valid_param_list'], 2] = d['l1_loss']
                val_loss_sets[val_dir_i * self.trained_epochs:(val_dir_i + 1) *
                              self.trained_epochs][d['valid_param_list'], 3] = d['ks_loss']
                valid_val_epoch_count += len(d['valid_param_list'])

            else:
                for epoch in range(self.trained_epochs):
                    param_val_path = os.path.join(val_dir, f'best_param{epoch}.pth')
                    if not os.path.exists(param_val_path):
                        continue
                    valid_val_epoch_count += 1
                    d = torch.load(param_val_path, map_location=torch.device(self.device))
                    parameter_sets[:, val_dir_i * self.trained_epochs + epoch] = torch.squeeze(d['parameter'])
                    val_loss_sets[val_dir_i * self.trained_epochs + epoch, 1] = d['corr_loss']
                    val_loss_sets[val_dir_i * self.trained_epochs + epoch, 2] = d['l1_loss']
                    val_loss_sets[val_dir_i * self.trained_epochs + epoch, 3] = d['ks_loss']
        if valid_val_dir_count == 0:
            print("No valid validated directories.")
            return 1
        if valid_val_epoch_count == 0:
            print("No valid epoch.")
            return 1
        val_loss_sets[:, 0] = torch.sum(val_loss_sets[:, 1:], dim=1)
        # Record all param_10 and loss
        val_total_loss, sorted_index = torch.sort(val_loss_sets[:, 0], descending=False)
        val_total_loss = val_total_loss[:self.param_sets]
        sorted_index = sorted_index[:self.param_sets]
        all_loss = val_loss_sets[sorted_index]
        parameter = parameter_sets[:, sorted_index]

        save_dict = {
            'parameter': parameter,
            'val_total_loss': val_total_loss,
            'sorted_index': sorted_index,
            'corr_loss': all_loss[:, 1],
            'l1_loss': all_loss[:, 2],
            'ks_loss': all_loss[:, 3]
        }
        torch.save(save_dict, os.path.join(self.test_dir, 'val_results.pth'))
        return 0

    def select_best_from_train(self):
        """Select the best performed parameter in training set, normally a step before E/I ratio calculation
        Particularly used for cases when there's only training set (like individual-level with only 1 run).

        Returns:
            int: state.
        """

        train_dirs = self.val_dirs
        train_dirs_len = len(train_dirs)
        parameter_sets = torch.zeros(self.parameters_dim, train_dirs_len * self.trained_epochs)  # [parameters_dim, 500]
        train_loss_sets = torch.ones(train_dirs_len * self.trained_epochs, 4) * 3
        # [total, corr, L1, ks]
        valid_train_dir_count = 0
        for train_dir_i in range(train_dirs_len):
            train_dir = train_dirs[train_dir_i]
            if not os.path.exists(f'{train_dir}/final_state_pFIC.pth'):
                print(f"{train_dir}/final_state.pth cannot be found.")
                continue
            valid_train_dir_count += 1
            for epoch in range(self.trained_epochs):
                param_save_path = os.path.join(train_dir, f'param_save_epoch{epoch}.pth')
                if not os.path.exists(param_save_path):
                    continue
                d = torch.load(param_save_path, map_location=torch.device(self.device))
                FC_FCD_loss = d['FC_FCD_loss']
                total_loss = torch.sum(FC_FCD_loss, dim=1)
                index_min_loss = torch.argmin(total_loss).item()
                index_min_parameter = d['valid_param_list'][index_min_loss]
                parameter_sets[:, train_dir_i * self.trained_epochs + epoch] = d[
                    'parameter'][:, index_min_parameter]  # d['parameter']:[parameter_dim, param_sets]
                train_loss_sets[train_dir_i * self.trained_epochs + epoch, 1:] = d['FC_FCD_loss'][index_min_loss]
        if valid_train_dir_count == 0:
            print("No valid directories.")
            return 1
        train_loss_sets[:, 0] = torch.sum(train_loss_sets[:, 1:], dim=1)
        # Record all param_10 and loss
        train_total_loss, sorted_index = torch.sort(train_loss_sets[:, 0], descending=False)
        train_total_loss = train_total_loss[:self.param_sets]
        sorted_index = sorted_index[:self.param_sets]
        sorted_loss = train_loss_sets[sorted_index]
        parameter = parameter_sets[:, sorted_index]

        save_dict = {
            'parameter': parameter,
            'train_total_loss': train_total_loss,
            'sorted_index': sorted_index,
            'corr_loss': sorted_loss[:, 1],
            'l1_loss': sorted_loss[:, 2],
            'ks_loss': sorted_loss[:, 3]
        }
        torch.save(save_dict, os.path.join(self.test_dir, 'train_results.pth'))
        return 0
