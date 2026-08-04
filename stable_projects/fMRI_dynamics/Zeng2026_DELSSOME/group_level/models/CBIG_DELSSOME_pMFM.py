"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import math
import torch
import numpy as np
import time
from tqdm import tqdm


def torch_corr_3D(vec_3d):
    """
    Compute the correlation coefficient parallel for 3D vector.
    And the 2nd dimension must be 2, standing for X and Y.
    :param vec_3d: [M, 2, len]
    :return: The corr_coef for X and Y. corr_coef: [M,]
    """
    ex_ey = torch.mean(vec_3d, dim=-1)
    std = torch.std(vec_3d, dim=-1, unbiased=False)
    xy = vec_3d[:, 0, :] * vec_3d[:, 1, :]
    e_xy = torch.mean(xy, dim=-1)
    cov = e_xy - ex_ey[:, 0] * ex_ey[:, 1]
    corr = cov / (std[:, 0] * std[:, 1])
    return corr


class MfmModel:
    """
    The MFM model class, which implements Deco 2013 and Xiaolu Kong 2021 model.
    """

    def __init__(self, config, parameter, sc_euler, dt):
        """
        :param parameter: (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    Each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent connection strength (w)
                    parameter[N:2*N]: external input current (I)
                    parameter[2*N]: Global constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        :param sc_euler: N*N structural connectivity matrix, should be normalized to max of 0.02
        :param dt: time interval
        """
        torch.set_default_dtype(torch.float64)

        self.N = (parameter.shape[0] - 1) // 3  # N = 68 ROIs
        self.M = parameter.shape[1]  # num of parameter sets
        self.parameter = parameter
        self.sc_euler = sc_euler
        self.dt = dt
        self.w = parameter[0:self.N]
        self.I_0 = parameter[self.N:2 * self.N]
        self.G = parameter[2 * self.N]
        self.sigma = parameter[2 * self.N + 1:3 * self.N + 1]

        # Synaptic Dynamical Equations Constants
        synaptic_constants = config['Dynamic Equation Constants']
        self.a = float(synaptic_constants['a'])
        self.b = float(synaptic_constants['b'])
        self.d = float(synaptic_constants['d'])
        self.tau_s = float(synaptic_constants['tau_s'])
        self.J = float(synaptic_constants['J'])
        self.gamma_kin = float(synaptic_constants['gamma_kin'])

        # Hemodynamic Model Constants
        hemodynamic_constants = config['Hemodynamic Model Constants']
        self.V0 = float(hemodynamic_constants['V0'])  # resting blood volume fraction
        self.kappa = float(hemodynamic_constants['kappa'])  # [s^-1] rate of signal decay
        self.gamma_hemo = float(hemodynamic_constants['gamma_hemo'])  # [s^-1] rate of flow-dependent elimination
        self.tau = float(hemodynamic_constants['tau'])  # [s] hemodynamic transit time
        self.alpha = float(hemodynamic_constants['alpha'])  # Grubb's exponent
        self.rho = float(hemodynamic_constants['rho'])  # resting oxygen extraction fraction
        self.B0 = float(hemodynamic_constants['B0'])  # magnetic field strength, T
        self.TE = float(hemodynamic_constants['TE'])  # TE echo time, s
        self.r0 = float(hemodynamic_constants['r0'])  # the intravascular relaxation rate, Hz
        self.epsilon_hemo = float(hemodynamic_constants['epsilon_hemo'])
        # the ratio between intravascular and extravascular MR signal
        self.k1 = 4.3 * 28.265 * self.B0 * self.TE * self.rho
        self.k2 = self.epsilon_hemo * self.r0 * self.TE * self.rho
        self.k3 = 1 - self.epsilon_hemo

        print("Successfully init MFM model!")

    def CBIG_mfm_simulation(self, simulate_time, burn_in_time, TR, warm_up_t, use_tqdm=False):
        """
        For M sets of parameters. Do not store the whole process like S_E, S_I to save memory.
        :param simulate_time: simulation time needed, in [min]
        :param burn_in_minute: The burn-in time (can be seen as the second warm up) in [min]
        :param TR: simulated time recording interval
        :param warm_up_t: The duration of warm up loop in frames
        :param use_tqdm: use progress bar to visualize simulation steps
        :param need_EI: return S_E_ave and S_I_ave
        :return: BOLD: [N, M, t]; valid_M_mask: [M], a boolean mask to indicate whether this parameter is valid
        """

        N = self.N
        M = self.M
        dt = self.dt

        # Set time
        t_start = 0  # seconds
        burn_in_t = 60 * burn_in_time
        t_end = t_start + burn_in_t + 60 * simulate_time
        t_p = torch.arange(t_start, t_end + dt, dt)
        t_len = len(t_p)
        t_inter = int(round(TR / dt))

        # Set initial values
        z = torch.zeros(N, M)
        f = torch.ones(N, M)
        v_volume = torch.ones(N, M)
        q = torch.ones(N, M)
        bold = torch.zeros(N, M, t_len // t_inter + 1)
        count_bold = 0
        S_t = torch.ones(N, M) * 0.001

        bold_time_start = time.time()
        print("Start BOLD calculating...")
        print("    Start warming up ...")
        if use_tqdm:
            warm_loop = tqdm(range(warm_up_t), position=0, leave=True)

            main_loop = tqdm(range(t_len), position=0, leave=True)
        else:
            warm_loop = range(warm_up_t)
            main_loop = range(t_len)

        for t1 in warm_loop:
            dS_dt = self._synaptic_dynamical_equations(S_t)
            S_t = S_t + dS_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)

        if use_tqdm:
            warm_loop.close()
        print("    End warming up and start main body...")

        S_t_ave = torch.zeros(N, M)

        # Start calculating
        for t in main_loop:
            dS_dt = self._synaptic_dynamical_equations(S_t)
            S_t = S_t + dS_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            # here math.sqrt(dt) is to make noise equivalent under different time interval, see notes.
            dz_dt, df_dt, dv_dt, dq_dt = self._hemodynamic_equations(S_t, z, f, v_volume, q)
            z = z + dz_dt * dt
            f = f + df_dt * dt
            v_volume = v_volume + dv_dt * dt
            q = q + dq_dt * dt

            S_t_ave = S_t_ave + S_t

            if (t + 2) % t_inter == 0:  # record every TR
                bold[:, :, count_bold] = 100 / self.rho * self.V0 * (self.k1 * (1 - q) + self.k2 *
                                                                     (1 - q / v_volume) + self.k3 * (1 - v_volume))
                # here 100 / rho is a scalar
                count_bold += 1

        bold[:, :, count_bold] = 100 / self.rho * self.V0 * (self.k1 * (1 - q) + self.k2 *
                                                             (1 - q / v_volume) + self.k3 * (1 - v_volume))
        # end main loop
        burn_in_bold = int(burn_in_t / TR)

        S_t_ave = S_t_ave / t_len

        # Check excitatory firing rate to filter those fMRI outside firing rate range.
        valid_M_mask = torch.ones(M, dtype=torch.bool)
        for i in range(self.M):
            if torch.isnan(bold[:, i, :]).any():
                valid_M_mask[i] = False

        bold_elapsed = time.time() - bold_time_start
        print('Time using for calculating BOLD signals cost: ', bold_elapsed)
        print("BOLD time series shape with burn-in: ", bold.shape)

        return bold[:, :, burn_in_bold + 1:], valid_M_mask

    def _synaptic_dynamical_equations(self, S_t):
        """
        The equations of Xiaolu 2021 model
        :param S_t: [N, M]
        :return: dS/dt
        """
        # Equations
        # Total input current
        I_t = self.w * self.J * S_t + self.G * self.J * torch.matmul(self.sc_euler, S_t) + self.I_0
        # firing rate
        r = (self.a * I_t - self.b) / (1 - torch.exp(-self.d * (self.a * I_t - self.b)))
        dS_dt = -S_t / self.tau_s + (1 - S_t) * self.gamma_kin * r
        return dS_dt

    def _hemodynamic_equations(self, S_t, z_t, f_t, v_t, q_t):
        dz_dt = S_t - self.kappa * z_t - self.gamma_hemo * (f_t - 1)
        df_dt = z_t
        dv_dt = (f_t - v_t**(1 / self.alpha)) / self.tau
        dq_dt = (f_t / self.rho * (1 - (1 - self.rho)**(1 / f_t)) - q_t * v_t**(1 / self.alpha - 1)) / self.tau
        return dz_dt, df_dt, dv_dt, dq_dt

    @staticmethod
    def all_loss_calculate_from_fc_fcd(fc_sim, fcd_hist_sim, fc_emp, fcd_cdf_emp):
        """
        Calculate corr_loss, L1_loss and KS_loss from simulated FC matrix and FCD histogram
        :param fc_sim: [M, N, N], sets of N x N simulated FC matrices
        :param fcd_hist: [fcd_hist_bins, M], no need to cumsum or normalize. Will do automatically in KS_cost
        :param fc_emp: [N, N], the empirical FC matrix
        :param fcd_cdf_emp: [fcd_hist_bins, 1], the empirical FCD CDF
        :return: [M], sets of loss
        """
        corr, L1_loss = MfmModel.FC_correlation_n_L1_cost(fc_sim, fc_emp)
        corr_loss = 1 - corr
        ks_loss = MfmModel.KS_cost(fcd_hist_sim, fcd_cdf_emp)
        total_loss = corr_loss + L1_loss + ks_loss
        return total_loss, corr_loss, L1_loss, ks_loss

    @staticmethod
    def all_loss_calculate_from_bold(bold, fc_emp, fcd_cdf_emp):
        """
        Calculate corr_loss, L1_loss and KS_loss from BOLD signals
        :param bold: [N, M, len], the simulated BOLD signals
        :param fc_emp: [N, N], the empirical FC matrix
        :param fcd_cdf_emp: [fcd_cdf_emp, 1]. Has been done cumulative summation and normalization
        (dividing by emp_fcd_cum[-1:, :])
        :return: total_loss: [M]
        """
        fc_sim = MfmModel.FC_calculate(bold)
        corr, L1_loss = MfmModel.FC_correlation_n_L1_cost(fc_sim, fc_emp)
        corr_loss = 1 - corr
        _, fcd_hist = MfmModel.FCD_calculate(bold)
        ks_loss = MfmModel.KS_cost(fcd_hist, fcd_cdf_emp)
        total_loss = corr_loss + L1_loss + ks_loss
        return total_loss

    @staticmethod
    def FC_calculate(bold):
        """
        Calculate FC matrix for M sets of BOLD signals
        :param bold: [N, M, t_len]
        :return: FC matrix [M, N, N]
        """
        N = bold.shape[0]
        M = bold.shape[1]
        fc_mat = torch.zeros((M, N, N))
        for i in range(M):
            fc_mat[i] = torch.corrcoef(bold[:, i, :])
        return fc_mat

    @staticmethod
    def FC_correlation_n_L1_cost(fc_sim, fc_emp):
        """
        Compute the FC correlation and L1 cost for all M sets.
        [ATTENTION] here L1 is not the real L1 (absolute difference of values and get the mean),
        but the absolute difference of two mean values.
        :param fc_sim: [M, N, N]
        :param fc_emp: [N, N]
        :return: corr: [M,]; L1: [M,]
        """
        M = fc_sim.shape[0]
        N = fc_sim.shape[1]

        # Extract the upper triangular part
        mask = torch.ones(N, N, dtype=torch.bool)
        mask = torch.triu(mask, 1)
        vec_emp = fc_emp[mask]
        vec_emp = vec_emp.unsqueeze(0).expand(M, -1)  # [M, len]
        vec_sim = fc_sim[:, mask]

        # L1 version 1: abs(mean)
        L1_cost = torch.abs(torch.mean(vec_emp, dim=1) - torch.mean(vec_sim, dim=1))

        vec_3d = torch.zeros(M, 2, vec_emp.shape[1])
        vec_3d[:, 0, :] = vec_sim
        vec_3d[:, 1, :] = vec_emp

        corr = torch_corr_3D(vec_3d)
        return corr, L1_cost

    @staticmethod
    def FC_correlation_single(fc_1, fc_2):
        """
        Compute the correlation between two FC matrix.
        By flatten their upright part to two vectors and use Pearson's correlation
        :param fc_1: [N, N]
        :param fc_2: [N, N]
        :return:
        """
        N = fc_1.shape[0]
        mask = np.ones((N, N)).astype(bool)
        mask = np.triu(mask, 1)
        vec = np.zeros((2, (N * N - N) // 2))
        vec[0] = fc_1[mask]
        vec[1] = fc_2[mask]
        cor_fc = np.corrcoef(vec)
        return cor_fc[0, 1]

    @staticmethod
    def FCD_calculate(bold, window_size=83, bins=10000):
        """
        Moving windows and calculating the FC matrix for every window_size,
        then calculating the correlation between these FC matrix.
        :param bold: [N, M, t_len]
        :param window_size: The length of sliding window
        :param bins: The histogram bins
        :return: FCD matrix [M, window_num, window_num]. window_num = t_len - window_size + 1
                FCD histogram: [bins, M]
        """
        if torch.cuda.is_available():
            torch.set_default_tensor_type('torch.cuda.DoubleTensor')
        else:
            torch.set_default_tensor_type('torch.DoubleTensor')

        N = bold.shape[0]
        M = bold.shape[1]
        t_len = bold.shape[2]
        window_num = t_len - window_size + 1
        if t_len < window_size:
            raise Exception("The length of bold signal is shorter than the window size!")

        # Compute sliding window FC
        fc_list = torch.zeros(M, window_num, N, N)
        for t in range(0, window_num):
            bold_single = bold[:, :, t:t + window_size]
            fc_list[:, t, :, :] = MfmModel.FC_calculate(bold_single)

        # Compute pairwise correlation between FCs
        fcd_mat = torch.zeros(M, window_num, window_num)
        fc_mask = torch.ones(N, N, dtype=torch.bool)
        fc_mask = torch.triu(fc_mask, 1)
        for m in range(M):
            fcd_mat[m] = torch.corrcoef(fc_list[m, :, fc_mask])

        # Calculate the FCD histogram
        fcd_mask = torch.ones(window_num, window_num, dtype=torch.bool)
        fcd_mask = torch.triu(fcd_mask, 1)
        fcd_vec = fcd_mat[:, fcd_mask]
        fcd_hist = torch.ones(bins, M)
        for hist_i in range(M):
            fcd_hist[:, hist_i] = torch.histc(fcd_vec[hist_i], bins=bins, min=-1., max=1.)

        return fcd_mat, fcd_hist

    @staticmethod
    def KS_cost(fcd_hist_sim, fcd_cdf_emp):
        """
        Calculate the KS cost between the simulated FCD histogram and the empirical one
        :param sim_fcd_hist: [bins, M]
        :param emp_fcd_cum: [bins, 1]. Has been done cumulative summation and normalization
        (divided by emp_fcd_cum[-1:, :])
        :return: KS_cost: [M]
        """
        M = fcd_hist_sim.shape[1]
        sim_fcd_cum = torch.cumsum(fcd_hist_sim, dim=0)
        sim_fcd_cum = sim_fcd_cum / sim_fcd_cum[-1:, :]
        emp_fcd_cum_expand = fcd_cdf_emp.expand(-1, M)
        ks_dif = torch.abs(sim_fcd_cum - emp_fcd_cum_expand)
        ks_cost = torch.max(ks_dif, dim=0)[0]
        return ks_cost
