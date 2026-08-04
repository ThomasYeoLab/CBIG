"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import math
import torch
import numpy as np
import time
from scipy.optimize import fsolve
from tqdm import tqdm


class MfmModel:

    def __init__(self, config, parameter, sc_euler, dt):
        """
        :param parameter: (N*3+1)*M matrix.
                    N is the number of ROI
                    M is the number of candidate parameter sets.
                    Each column of matrix presents a parameter set, where:
                    parameter[0:N]: recurrent strength within excitatory population (wEE)
                    parameter[N:2*N]: connection strength from excitatory population to inhibitory population (wEI)
                    parameter[2*N]: Global constant G
                    parameter[2*N+1:3*N+1]: noise amplitude sigma
        :param sc_euler: N*N structural connectivity matrix, should be normalized to max of 0.02
        :param dt: time interval
        """

        torch.set_default_dtype(torch.float64)

        # Model Adaptation
        self.is_w_IE_fixed = True

        self.N = (parameter.shape[0] - 1) // 3  # N = 68 ROIs
        self.M = parameter.shape[1]  # num of parameter sets
        self.parameter = parameter
        self.sc_euler = sc_euler
        self.dt = dt
        self.w_EE = parameter[0:self.N]
        self.w_EI = parameter[self.N:2 * self.N]
        self.G = parameter[2 * self.N]
        self.sigma = parameter[2 * self.N + 1:3 * self.N + 1]

        # Synaptic Dynamical Equations Constants
        synaptic_constants = config['Dynamic Equation Constants']
        self.I_0 = float(synaptic_constants['I_0'])
        self.a_E = float(synaptic_constants['a_E'])
        self.b_E = float(synaptic_constants['b_E'])
        self.d_E = float(synaptic_constants['d_E'])
        self.tau_E = float(synaptic_constants['tau_E'])
        self.W_E = float(synaptic_constants['W_E'])
        self.a_I = float(synaptic_constants['a_I'])
        self.b_I = float(synaptic_constants['b_I'])
        self.d_I = float(synaptic_constants['d_I'])
        self.tau_I = float(synaptic_constants['tau_I'])
        self.W_I = float(synaptic_constants['W_I'])
        self.J_NMDA = float(synaptic_constants['J_NMDA'])
        self.gamma_kin = float(synaptic_constants['gamma_kin'])
        self.r_E = float(synaptic_constants['r_E'])

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

        # Solve for steady-state synaptic variables (S_E, I_E, I_I, w_IE) using closed-form fixed-point equations
        # from the MFM2013 parameterization before starting the simulation.
        # Calculate S_E average
        S_E_ini = 0.1641205151
        self.S_E_ave, info_dict, ier, message = fsolve(self._solve_S_E_ave, S_E_ini, full_output=True)
        if ier == 0:
            print(message)
        self.S_E_ave = self.S_E_ave[0]  # convert 1 element array to float

        # Calculate I_E average
        I_E_ini = 0.3772259651
        I_E_ave, info_dict, ier, message = fsolve(self._solve_I_E_ave, I_E_ini, full_output=True)
        if ier == 0:
            print(message)
        I_E_ave = I_E_ave[0]

        # Calculate I_I fixed point
        I_I_ini = 0.296385800197336 * np.ones((self.N, self.M))
        I_I_ave = np.ones_like(I_I_ini)
        for m in range(self.M):
            I_I_func_args = self.S_E_ave, self.w_EI[:, m]
            I_I_ave[:, m], info_dict, ier, message = fsolve(self._solve_I_I,
                                                            I_I_ini[:, m],
                                                            args=I_I_func_args,
                                                            full_output=True)
            if ier == 0:
                print(message)

        # Calculate w_IE
        S_I_ave = self.tau_I * (self.a_I * I_I_ave - self.b_I) / (1 - np.exp(-self.d_I *
                                                                             (self.a_I * I_I_ave - self.b_I)))
        S_I_ave = torch.as_tensor(S_I_ave)
        self.S_I_ave = S_I_ave
        self.w_IE = (self.W_E * self.I_0 + self.w_EE * self.J_NMDA * self.S_E_ave + self.G * self.J_NMDA *
                     torch.matmul(self.sc_euler, self.S_E_ave * torch.ones(self.N, self.M)) - I_E_ave) / S_I_ave
        # w_IE: [N, M]

        print("pFIC model successfully initialized... ")

    def CBIG_mfm_simulation(self, simulate_time, burn_in_time, TR, warm_up_t, use_tqdm=False, need_EI=False):
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
        r_E = torch.zeros(N, M)
        z = torch.zeros(N, M)
        f = torch.ones(N, M)
        v_volume = torch.ones(N, M)
        q = torch.ones(N, M)
        bold = torch.zeros(N, M, t_len // t_inter + 1)
        count_bold = 0
        S_E = torch.ones(N, M) * self.S_E_ave
        S_I = torch.ones(N, M) * 0.1433408985

        bold_time_start = time.time()
        print("Start warm-up ...")
        if use_tqdm:
            warm_loop = tqdm(range(warm_up_t), position=0, leave=True)

            main_loop = tqdm(range(t_len), position=0, leave=True)
        else:
            warm_loop = range(warm_up_t)
            main_loop = range(t_len)

        for t1 in warm_loop:
            dSE_dt, dSI_dt, _ = self._synaptic_dynamical_equations(S_E, S_I)
            S_E = S_E + dSE_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            S_I = S_I + dSI_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)

        if use_tqdm:
            warm_loop.close()
        print("Start forward simulation...")

        S_E_ave = torch.zeros(N, M)
        S_I_ave = torch.zeros(N, M)
        r_E_ave = torch.zeros(N, M)

        # Start calculating
        for t in main_loop:
            dSE_dt, dSI_dt, r_E = self._synaptic_dynamical_equations(S_E, S_I)
            S_E = S_E + dSE_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            # here math.sqrt(dt) is to make noise equivalent under different time interval, see notes.
            S_I = S_I + dSI_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            dz_dt, df_dt, dv_dt, dq_dt = self._hemodynamic_equations(S_E, z, f, v_volume, q)
            z = z + dz_dt * dt
            f = f + df_dt * dt
            v_volume = v_volume + dv_dt * dt
            q = q + dq_dt * dt

            S_E_ave = S_E_ave + S_E
            S_I_ave = S_I_ave + S_I
            r_E_ave = r_E_ave + r_E

            if (t + 2) % t_inter == 0:  # record every TR
                bold[:, :, count_bold] = 100 / self.rho * self.V0 * (self.k1 * (1 - q) + self.k2 *
                                                                     (1 - q / v_volume) + self.k3 * (1 - v_volume))
                # here 100 / rho is a scalar
                count_bold += 1

        bold[:, :, count_bold] = 100 / self.rho * self.V0 * (self.k1 * (1 - q) + self.k2 *
                                                             (1 - q / v_volume) + self.k3 * (1 - v_volume))
        # end main loop
        burn_in_bold = int(burn_in_t / TR)

        S_E_ave = S_E_ave / t_len
        S_I_ave = S_I_ave / t_len
        r_E_ave = r_E_ave / t_len

        # Filter out parameter sets where the simulation diverged (NaN in firing rate or BOLD signal);
        # downstream loss is only computed on valid runs.
        valid_M_mask = torch.ones(M, dtype=torch.bool)
        for i in range(self.M):
            if torch.isnan(r_E_ave[:, i]).any():
                valid_M_mask[i] = False
            elif torch.isnan(bold[:, i, :]).any():
                valid_M_mask[i] = False

        bold_elapsed = time.time() - bold_time_start
        print('Time eplased: ', bold_elapsed)

        if need_EI:
            return bold[:, :, burn_in_bold + 1:], valid_M_mask, S_E_ave, S_I_ave
        else:
            return bold[:, :, burn_in_bold + 1:], valid_M_mask

    def _synaptic_dynamical_equations(self, S_E_t, S_I_t):
        """
        The equations of Deco 2014 model
        :param S_E_t: [N, M]
        :param S_I_t: [N, M]
        :return: dS^E/dt, dS^I/dt
        """
        # Equations
        '''
        if self.is_w_IE_fixed:
            w_IE = self.w_IE
        else:
            w_IE = 0
        '''
        I_E_t = self.W_E * self.I_0 + self.w_EE * self.J_NMDA * S_E_t + self.G * self.J_NMDA * \
            torch.matmul(self.sc_euler, S_E_t) - self.w_IE * S_I_t
        I_I_t = self.W_I * self.I_0 + self.w_EI * self.J_NMDA * S_E_t - S_I_t
        r_E = (self.a_E * I_E_t - self.b_E) / (1 - torch.exp(-self.d_E * (self.a_E * I_E_t - self.b_E)))
        r_I = (self.a_I * I_I_t - self.b_I) / (1 - torch.exp(-self.d_I * (self.a_I * I_I_t - self.b_I)))
        dSE_dt = -S_E_t / self.tau_E + (1 - S_E_t) * self.gamma_kin * r_E
        dSI_dt = -S_I_t / self.tau_I + r_I
        return dSE_dt, dSI_dt, r_E

    def _hemodynamic_equations(self, S_E_t, z_t, f_t, v_t, q_t):
        dz_dt = S_E_t - self.kappa * z_t - self.gamma_hemo * (f_t - 1)
        df_dt = z_t
        dv_dt = (f_t - v_t**(1 / self.alpha)) / self.tau
        dq_dt = (f_t / self.rho * (1 - (1 - self.rho)**(1 / f_t)) - q_t * v_t**(1 / self.alpha - 1)) / self.tau
        return dz_dt, df_dt, dv_dt, dq_dt

    def _solve_I_I(self, I_I_ave, S_E_ave, w_EI_m):
        if torch.cuda.is_available():
            w_EI_m = w_EI_m.cpu().numpy()
        phi_I_I_ave = (self.a_I * I_I_ave - self.b_I) / (1 - np.exp(-self.d_I * (self.a_I * I_I_ave - self.b_I)))
        res = self.W_I * self.I_0 + self.J_NMDA * w_EI_m * S_E_ave - phi_I_I_ave * self.tau_I - I_I_ave
        return res

    def _solve_S_E_ave(self, S_E_ave):
        res = S_E_ave / (self.tau_E * self.gamma_kin * (1 - S_E_ave)) - self.r_E
        # Fixed-point equation for S_E_ave: S_E_ave = 0.0641 r_E / (1 + 0.0641 r_E).
        # Kept in this exact form to reproduce the published results.
        return res

    def _solve_I_E_ave(self, I_E_ave):
        tmp = self.a_E * I_E_ave - self.b_E
        res = tmp / (1 - np.exp(-self.d_E * tmp)) - self.r_E
        return res
