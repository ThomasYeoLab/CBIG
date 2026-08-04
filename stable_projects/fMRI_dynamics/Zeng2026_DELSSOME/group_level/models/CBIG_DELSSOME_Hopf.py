"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
import torch
import time
from tqdm import tqdm
import math


class HopfModel:
    """
    Class for Hopf model.

    Attributes
    ----------
    params : dict
        Parameters of the Hopf model.
    """

    def __init__(self, config, parameter, sc_euler, dt):
        """
        Initialize the Hopf model with given parameters.

        Parameters
        ----------
        params : dict
            Parameters for the Hopf model.
        """
        torch.set_default_dtype(torch.float64)
        self.N = (parameter.shape[0] - 1) // 3  # N = 68 ROIs
        self.M = parameter.shape[1]  # num of parameter sets
        self.parameter = parameter
        self.sc_euler = sc_euler
        self.dt = dt
        self.a = parameter[0:self.N]  # [N, M]
        self.omega = parameter[self.N:2 * self.N]
        self.G = parameter[2 * self.N]
        self.sigma = parameter[2 * self.N + 1:3 * self.N + 1]

        print("Successfully init Hopf model.")

    def simulate(self, simulation_period, burn_in_period, TR, warmup_period, use_tqdm=False):
        """
        Simulate the Hopf model.

        Parameters
        ----------
        initial_conditions : list
            Initial conditions for the simulation.
        time_span : tuple
            Time span for the simulation.

        Returns
        -------
        results : list
            Simulation results.
        """
        N = self.N
        M = self.M
        dt = self.dt

        # Set time
        t_start = 0  # seconds
        burn_in_t = 60 * burn_in_period  # seconds
        t_end = t_start + burn_in_t + 60 * simulation_period  # seconds
        t_p = torch.arange(t_start, t_end + dt, dt)
        t_len = len(t_p)
        t_inter = int(round(TR / dt))

        # Set initial values
        x_real_ts = torch.zeros(N, M, t_len // t_inter + 1)
        y_imag_ts = torch.zeros(N, M, t_len // t_inter + 1)
        x_real = torch.zeros(N, M)
        y_imag = torch.zeros(N, M)

        count_tp = 0

        start_time = time.time()
        print("Start time series calculating...")
        print("    Start warming up ...")
        if use_tqdm:
            warm_loop = tqdm(range(warmup_period), position=0, leave=True)

            main_loop = tqdm(range(t_len), position=0, leave=True)
        else:
            warm_loop = range(warmup_period)
            main_loop = range(t_len)

        for t1 in warm_loop:
            dx_dt, dy_dt = self._dynamical_equations(x_real, y_imag)
            x_real = x_real + dx_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            y_imag = y_imag + dy_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)

        if use_tqdm:
            warm_loop.close()
        print("    End warming up and start main body...")

        # Start calculating
        for t in main_loop:
            dx_dt, dy_dt = self._dynamical_equations(x_real, y_imag)
            x_real = x_real + dx_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            y_imag = y_imag + dy_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)

            if (t + 1) % t_inter == 0:  # record every TR
                x_real_ts[:, :, count_tp] = x_real
                y_imag_ts[:, :, count_tp] = y_imag
                count_tp += 1

        if use_tqdm:
            main_loop.close()
        # end main loop
        burn_in_ts = int(burn_in_t / TR)

        # Check excitatory firing rate to filter those fMRI outside firing rate range.
        valid_M_mask = torch.ones(M, dtype=torch.bool)
        for i in range(self.M):
            if torch.isnan(x_real_ts[:, i, :]).any():
                valid_M_mask[i] = False

        elapsed_time = time.time() - start_time
        print('Time using for calculating time series cost: ', elapsed_time)
        print("Time series shape with burn-in: ", x_real_ts.shape)

        return x_real_ts[:, :, burn_in_ts + 1:], valid_M_mask

    def _dynamical_equations(self, x_real, y_imag):
        """
        Calculate the dynamical equations for the Hopf model.

        Parameters
        ----------
        x_real : torch.Tensor
            Real part of the state variable.
        y_imag : torch.Tensor
            Imaginary part of the state variable.
        Returns
        -------
        dx_dt : torch.Tensor
            Derivative of the real part of the state variable.
        dy_dt : torch.Tensor
            Derivative of the imaginary part of the state variable.
        """
        dx_dt = (self.a - torch.square(x_real) - torch.square(y_imag)) * x_real - self.omega * y_imag + self.G * (
            self.sc_euler @ x_real - torch.sum(self.sc_euler, dim=1, keepdim=True) * x_real)
        dy_dt = (self.a - torch.square(x_real) - torch.square(y_imag)) * y_imag + self.omega * x_real + self.G * (
            self.sc_euler @ y_imag - torch.sum(self.sc_euler, dim=1, keepdim=True) * y_imag)
        return dx_dt, dy_dt
