# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Feedback-inhibition-control (FIC) mean-field model + Euler integration.

Implements the Deco et al., 2014 model plus the hemodynamic forward model.
Scope of this release:

* Only ``pFIC`` is included (MFM-2013 and the legacy within-range classifiers
  are not part of this release; the individual-level paper results all use
  pFIC).
* :meth:`pFIC.load_dl_models` no longer hard-codes cluster paths — instead
  the checkpoint path is provided via ``emp_stats['delssome_ckpt']`` so users
  can plug in their own trained DELSSOME model.
"""

import math
import datetime
import time

import torch
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm

from DELSSOME_indiv.models.CBIG_DELSSOME_predictor import MulModLossPredictorMfm
from DELSSOME_indiv.CBIG_misc_utils import set_torch_default, CBIG_corr_torch
from DELSSOME_indiv.CBIG_FC_FCD_utils import calc_FC, calc_FCD, calc_all_loss_from_fc_fcd


# --------------------------------------------------------------------------- #
# Reshape helpers
# --------------------------------------------------------------------------- #
def reshape_dup_sim_res(sim_res_dict, param_sets, param_dup, exclude_keys=None):
    """Separate the ``param_sets * param_dup`` axis into two named axes."""
    exclude_keys = exclude_keys or []
    reshaped = {}
    for key, value in sim_res_dict.items():
        if key in exclude_keys or not isinstance(value, torch.Tensor):
            reshaped[key] = value
            continue
        if value.ndim == 1:
            reshaped[key] = value.view(param_dup, param_sets).T
        elif value.ndim == 2:
            N = value.shape[0]
            reshaped[key] = value.view(N, param_dup, param_sets).transpose(1, 2)
        elif value.ndim == 3:
            N = value.shape[0]
            T = value.shape[2]
            reshaped[key] = value.view(N, param_dup, param_sets, T).transpose(1, 2)
        else:
            print(f"Warning: Key '{key}' has unsupported shape {value.shape}; "
                  "keeping it unchanged.")
            reshaped[key] = value
    return reshaped


def get_valid_params_after_dup(valid_M_mask):
    """Param-set is valid if at least one of its duplicates is valid."""
    if valid_M_mask.ndim == 1:
        valid_M_mask = valid_M_mask.unsqueeze(0)
    valid_param_indices = []
    for i in range(valid_M_mask.shape[0]):
        if valid_M_mask[i].any():
            valid_param_indices.append(i)
    return valid_param_indices


# --------------------------------------------------------------------------- #
# Base class
# --------------------------------------------------------------------------- #
class DynamicModel:
    """Abstract base for biophysical dynamic models."""

    def __init__(self, config, parameter, sc_euler, dt):
        self.device, self.tensor_type = set_torch_default()
        self.config = config
        self.parameter = parameter
        self.parameter_dim = parameter.shape[0]  # (N*3+1)
        self.M = parameter.shape[1]
        self.sc_euler = sc_euler
        self.dt = dt
        self._init_from_config()

    def _init_from_config(self):
        raise NotImplementedError

    @classmethod
    def get_parameter_dim(cls, N):
        return N * 3 + 1

    @property
    def N(self):
        return (self.parameter.shape[0] - 1) // 3

    # The interface below is implemented by pFIC (subclasses can override).
    @classmethod
    def get_search_range(cls, N, stats):
        raise NotImplementedError

    @classmethod
    def get_parameters(cls, param_coef, concat_mat):
        raise NotImplementedError

    @classmethod
    def is_parameter_valid(cls, config, emp_stats, parameter, search_range):
        raise NotImplementedError

    @classmethod
    def get_parameterization(cls, emp_stats, other_parameterization=None):
        raise NotImplementedError

    @classmethod
    def init_CMA_ES(cls, search_range, pinv_concat_mat, param_coef_dim):
        raise NotImplementedError

    @classmethod
    def load_dl_models(cls, config, emp_stats, device):
        raise NotImplementedError

    @classmethod
    def sim_n_get_losses_dl(cls, dl_models, param_vectors, emp_stats,
                            sort_param_vectors=False):
        raise NotImplementedError

    def simulate(self, simulate_time, burn_in_time, TR, warm_up_t=5000,
                 use_tqdm=False, **kwargs) -> dict:
        raise NotImplementedError

    def derive_more_dup_sim_res(self, dup_sim_res, **kwargs):
        print("No more derived simulation results for this model.")
        return dup_sim_res, None

    @classmethod
    def calc_losses(cls, save_sim_res, other_sim_res, emp_stats, apply_sort=True):
        raise NotImplementedError


# --------------------------------------------------------------------------- #
# pFIC — feedback inhibition control model (Deco et al., 2014)
# --------------------------------------------------------------------------- #
class pFIC(DynamicModel):
    """
    Implementation of the FIC mean-field model of Deco et al. (2014).

    The parameter vector has 3N+1 entries packed as ``[wEE; wEI; G; sigma]``.
    The model integrates SE/SI synaptic gating + Balloon-Windkessel hemodynamics
    using the Euler-Maruyama scheme with step ``dt``.
    """

    def _init_from_config(self):
        N = self.N
        self.is_w_IE_fixed = True

        # Unpack parameter vector
        self.w_EE = self.parameter[0:N]
        self.w_EI = self.parameter[N:2 * N]
        self.G = self.parameter[2 * N]
        self.sigma = self.parameter[2 * N + 1:3 * N + 1]

        # Biophysical constants (read from the .ini config)
        dc = self.config["Dynamic Model Constants"]
        self.I_0 = float(dc["I_0"])
        self.a_E = float(dc["a_E"])
        self.b_E = float(dc["b_E"])
        self.d_E = float(dc["d_E"])
        self.tau_E = float(dc["tau_E"])
        self.W_E = float(dc["W_E"])
        self.a_I = float(dc["a_I"])
        self.b_I = float(dc["b_I"])
        self.d_I = float(dc["d_I"])
        self.tau_I = float(dc["tau_I"])
        self.W_I = float(dc["W_I"])
        self.J_NMDA = float(dc["J_NMDA"])
        self.gamma_kin = float(dc["gamma_kin"])
        self.r_E = float(dc["r_E"])
        print(f"Excitatory firing rate will be around {self.r_E} Hz.")
        self.r_E_min = self.r_E - 0.3
        self.r_E_max = self.r_E + 0.3

        # Hemodynamic / BOLD constants
        hc = self.config["Hemodynamic Model Constants"]
        self.V0 = float(hc["V0"])
        self.kappa = float(hc["kappa"])
        self.gamma_hemo = float(hc["gamma_hemo"])
        self.tau = float(hc["tau"])
        self.alpha = float(hc["alpha"])
        self.rho = float(hc["rho"])
        self.B0 = float(hc["B0"])
        self.TE = float(hc["TE"])
        self.r0 = float(hc["r0"])
        self.epsilon_hemo = float(hc["epsilon_hemo"])
        self.k1 = 4.3 * 28.265 * self.B0 * self.TE * self.rho
        self.k2 = self.epsilon_hemo * self.r0 * self.TE * self.rho
        self.k3 = 1 - self.epsilon_hemo

        self.fcd_hist_bins = int(self.config["Simulating Parameters"]["fcd_hist_bins"])
        self.window_size = int(self.config["Dataset Parameters"]["window_size"])

        # Solve for fixed-point inhibitory current, then derive w_IE so that the
        # FIC condition (Deco et al. 2014, Sec. 2.2) is satisfied identically.
        S_E_ini = 0.1641205151
        self.S_E_ave, _, ier, msg = fsolve(self._solve_S_E_ave, S_E_ini,
                                           full_output=True)
        if ier == 0:
            print(msg)
        self.S_E_ave = self.S_E_ave[0]

        I_E_ini = 0.3772259651
        I_E_ave, _, ier, msg = fsolve(self._solve_I_E_ave, I_E_ini,
                                      full_output=True)
        if ier == 0:
            print(msg)
        I_E_ave = I_E_ave[0]

        # Per-parameter-set fixed-point I_I and the resulting S_I_ave.
        I_I_ini = 0.296385800197336 * np.ones((N, self.M))
        I_I_ave = np.ones_like(I_I_ini)
        for m in range(self.M):
            args = (self.S_E_ave, self.w_EI[:, m])
            I_I_ave[:, m], _, ier, msg = fsolve(
                self._solve_I_I, I_I_ini[:, m], args=args, full_output=True
            )
            if ier == 0:
                print(msg)

        S_I_ave = self.tau_I * (self.a_I * I_I_ave - self.b_I) / (
            1 - np.exp(-self.d_I * (self.a_I * I_I_ave - self.b_I))
        )
        S_I_ave = torch.as_tensor(S_I_ave).to(self.tensor_type)

        self.w_IE = (
            self.W_E * self.I_0
            + self.w_EE * self.J_NMDA * self.S_E_ave
            + self.G * self.J_NMDA * torch.matmul(
                self.sc_euler, self.S_E_ave * torch.ones(N, self.M))
            - I_E_ave
        ) / S_I_ave
        print("Successfully init MFM 2014 model!")

    # ---------------------------------------------------------------- #
    # Static methods used by CMA-ES (search range, parameterisation, ...)
    # ---------------------------------------------------------------- #
    @classmethod
    def get_search_range(cls, N, emp_stats):
        search_range = torch.zeros(cls.get_parameter_dim(N), 2)
        search_range[0:N, 0] = 1
        search_range[0:N, 1] = 10
        wei_min, wei_max = cls.get_wei_range(N, emp_stats)
        search_range[N:2 * N, 0] = wei_min
        search_range[N:2 * N, 1] = wei_max
        search_range[2 * N, 0] = 0
        search_range[2 * N, 1] = 3
        search_range[2 * N + 1:, 0] = 0.0005
        search_range[2 * N + 1:, 1] = 0.01
        return search_range

    @classmethod
    def get_wei_range(cls, N, emp_stats):
        """Heuristic wEI search range; tightened by empirical FC mean.

        At individual-level the FC tends to be noisier, so we cap wEI_max at 5
        and clip wEI_min to ``[0, 3.5]``, as used for the paper's results.
        """
        if emp_stats["group_type"] is not None and emp_stats["ds_name"] == "HCP-YA":
            wei_min, wei_max = 0, 5
        else:
            upper_trig_mask = torch.triu(torch.ones(N, N, dtype=torch.bool), 1)
            mean_fc = torch.mean(emp_stats["emp_fc"][upper_trig_mask])
            wei_max = 11 * mean_fc
            wei_max = min(max(wei_max, 3.2), 5)
            wei_min = 2 * wei_max - 6
            wei_min = min(max(wei_min, 0), 3.5)
        print(f"wEI search range [{wei_min}, {wei_max}].")
        return wei_min, wei_max

    @classmethod
    def get_parameters(cls, param_coef, concat_mat):
        """``param_coef -> 3N+1 FIC parameters`` via the myelin/RSFC concat matrix."""
        p = concat_mat.shape[1]
        w_EE = torch.matmul(concat_mat, param_coef[0:p])
        w_EI = torch.matmul(concat_mat, param_coef[p:2 * p])
        G = param_coef[2 * p].unsqueeze(0)
        sigma = torch.matmul(concat_mat, param_coef[2 * p + 1:])
        return torch.cat((w_EE, w_EI, G, sigma), dim=0).squeeze()

    @classmethod
    def is_parameter_valid(cls, config, emp_stats, parameter, search_range):
        sim_params = config["Simulating Parameters"]
        N = int(sim_params["n_ROI"])
        wEI_threshold = sim_params.get("wEI_variability_threshold", None)
        if wEI_threshold is not None:
            wEI_threshold = float(wEI_threshold)
        if wEI_threshold == 0:
            wEI_threshold = None

        outside_range = ((parameter < search_range[:, 0]).any()
                         or (parameter > search_range[:, 1]).any())
        if outside_range:
            return False
        if wEI_threshold is not None:
            wEI = parameter[N:2 * N].unsqueeze(1)
            wEI_myelin_corr = torch.squeeze(CBIG_corr_torch(wEI, emp_stats["myelin"]))
            wEI_rsfc_corr = torch.squeeze(CBIG_corr_torch(wEI, emp_stats["rsfc_gradient"]))
            # wEI should be POSITIVELY correlated with the RSFC gradient and
            # NEGATIVELY correlated with myelin (sensory-association axis).
            violates = (wEI_myelin_corr > 0) or (wEI_rsfc_corr < 0)
            too_uniform = torch.std(wEI) < wEI_threshold
            if violates or too_uniform:
                return False
        return True

    @classmethod
    def get_parameterization(cls, emp_stats, other_parameterization=None):
        """Build the concat matrix ``[1, myelin, rsfc_gradient]`` used to expand
        the 10 CMA-ES parameter-coefficients into N-vectors of ``wEE``, ``wEI``, ``sigma``."""
        myelin = emp_stats["myelin"]
        rsfc_gradient = emp_stats["rsfc_gradient"]
        if other_parameterization is None:
            concat_mat = torch.hstack(
                (torch.ones_like(myelin), myelin, rsfc_gradient))
            p = 3
        else:
            concat_mat = other_parameterization
            p = other_parameterization.shape[1]
        pinv_concat_mat = torch.linalg.pinv(concat_mat)
        param_coef_dim = 3 * p + 1  # wEE coeffs + wEI coeffs + G + sigma coeffs
        return param_coef_dim, concat_mat, pinv_concat_mat

    @classmethod
    def init_CMA_ES(cls, search_range, pinv_concat_mat, param_coef_dim):
        p, N = pinv_concat_mat.shape
        parameter_dim = cls.get_parameter_dim(N)

        # Initial parameter draw inside the search range.
        init_parameters = (torch.rand(parameter_dim)
                           * (search_range[:, 1] - search_range[:, 0])
                           + search_range[:, 0])
        init_parameters = init_parameters.unsqueeze(1)

        # Project per-region inits down to the 10-D coefficient space.
        start_point_wEE = torch.matmul(pinv_concat_mat, init_parameters[0:N]).squeeze()
        start_point_wEI = torch.matmul(pinv_concat_mat, init_parameters[N:2 * N]).squeeze()
        start_point_sigma = torch.matmul(pinv_concat_mat, init_parameters[2 * N + 1:]).squeeze()

        m_0 = torch.zeros(param_coef_dim)
        m_0[0:p] = start_point_wEE
        m_0[p:2 * p] = start_point_wEI
        m_0[2 * p] = init_parameters[2 * N]      # G
        m_0[2 * p + 1:] = start_point_sigma

        sigma_0 = 0.2
        p_sigma_0 = torch.zeros(param_coef_dim, 1)
        p_c_0 = torch.zeros(param_coef_dim, 1)
        V_ini = torch.eye(param_coef_dim)
        Lambda_ini = torch.ones(param_coef_dim)
        Lambda_ini[0:p] = start_point_wEE[0]
        Lambda_ini[p:2 * p] = start_point_wEI[0]
        Lambda_ini[2 * p] = 0.4
        Lambda_ini[2 * p + 1:] = 0.0005
        cov_0 = torch.matmul(V_ini,
                             torch.matmul(torch.diag(Lambda_ini ** 2), V_ini.T))
        return m_0, sigma_0, cov_0, p_sigma_0, p_c_0

    @classmethod
    def load_dl_models(cls, config, emp_stats, device):
        """
        Load the DELSSOME cost predictor used by DELSSOME CMA-ES.

        We expect ``emp_stats['delssome_ckpt']`` to point at a Lightning
        checkpoint of :class:`MulModLossPredictorMfm`. The within-range
        classifier is intentionally NOT used (Section 4.2 of the paper).
        """
        ckpt_path = emp_stats.get("delssome_ckpt", None)
        if ckpt_path is None:
            raise ValueError(
                "emp_stats['delssome_ckpt'] must point at a trained DELSSOME "
                "checkpoint (a Lightning .ckpt file from "
                "CBIG_DELSSOME_indiv_predictor_trainer.py)."
            )
        predictor = MulModLossPredictorMfm.load_from_checkpoint(
            ckpt_path, map_location=device,
        )
        predictor.eval()
        print("Using predictor:", predictor.__class__.__name__)
        dl_models = {"predictor": predictor, "is_multimodal": True}
        print("Successfully loaded DELSSOME predictor.")
        return dl_models

    @classmethod
    def sim_n_get_losses_dl(cls, dl_models, param_vectors, emp_stats,
                            sort_param_vectors=False):
        """
        Predict FC+FCD costs for a batch of candidate parameter vectors using
        DELSSOME instead of Euler integration. Returns the standard save-dict
        with ``total_loss``, ``corr_loss``, ``l1_loss``, ``ks_loss``, etc.
        """
        predictor = dl_models["predictor"]
        is_multimodal = dl_models.get("is_multimodal", True)
        classifier = dl_models.get("classifier", None)

        print("Simulating using deep learning models...")
        N = (param_vectors.shape[0] - 1) // 3
        with torch.no_grad():
            # ``param_vectors`` was [3N+1, M]; the DL model expects [M, 3N+1].
            param_vectors = param_vectors.T
            valid_param_indices = torch.arange(0, param_vectors.shape[0], 1)
            if classifier is not None:
                sc_for_a_param = emp_stats["sc_dl"]
                sc_for_all = sc_for_a_param.unsqueeze(0).expand(
                    param_vectors.shape[0], *sc_for_a_param.shape
                )
                pred_nan = torch.squeeze(classifier(param_vectors, sc_for_all))
                valid_param_indices = valid_param_indices[pred_nan > 0.5]
            valid_parameter = param_vectors[valid_param_indices]

            sc_for_a_param = emp_stats["sc_euler" if is_multimodal else "sc_dl"]
            sc_for_valid = sc_for_a_param.unsqueeze(0).expand(
                valid_parameter.shape[0], *sc_for_a_param.shape)
            fc_for_a_param = emp_stats["emp_fc" if is_multimodal else "fc_dl"]
            fc_for_valid = fc_for_a_param.unsqueeze(0).expand(
                valid_parameter.shape[0], *fc_for_a_param.shape)
            fcd_for_a_param = emp_stats["fcd_dl"]
            fcd_for_valid = fcd_for_a_param.unsqueeze(0).expand(
                valid_parameter.shape[0], *fcd_for_a_param.shape)

            if is_multimodal:
                G_valid = valid_parameter[:, 2 * N]
                # Scale every region's SC by the candidate G value.
                sc_for_valid = sc_for_valid * G_valid.unsqueeze(-1).unsqueeze(-1)
                # Reshape parameter into per-ROI tokens [M, R, 3].
                valid_parameter = torch.stack(
                    [
                        valid_parameter[:, :N],
                        valid_parameter[:, N:2 * N],
                        valid_parameter[:, 2 * N + 1:],
                    ],
                    dim=2,
                )

            pred_loss = predictor(valid_parameter, sc_for_valid,
                                  fc_for_valid, fcd_for_valid)  # [M_valid, 3]
            total_loss = torch.sum(pred_loss, dim=1)
            corr_loss = pred_loss[:, 0]
            l1_loss = pred_loss[:, 1]
            ks_loss = pred_loss[:, 2]

            # Sort small -> large total loss.
            total_loss, idx_sorted = torch.sort(total_loss, descending=False)
            valid_parameter = valid_parameter[idx_sorted]
            valid_param_indices = valid_param_indices[idx_sorted]
            corr_loss = corr_loss[idx_sorted]
            l1_loss = l1_loss[idx_sorted]
            ks_loss = ks_loss[idx_sorted]
            if sort_param_vectors:
                param_vectors = param_vectors[valid_param_indices]
            return {
                "total_loss": total_loss,
                "corr_loss": corr_loss,
                "l1_loss": l1_loss,
                "ks_loss": ks_loss,
                "parameter": param_vectors.T,
                "valid_param_indices": valid_param_indices,
            }

    # ---------------------------------------------------------------- #
    # Euler integration
    # ---------------------------------------------------------------- #
    def simulate(self, simulate_time, burn_in_time, TR, warm_up_t=5000,
                 use_tqdm=False, **kwargs) -> dict:
        """
        Euler-Maruyama integration of SE/SI + Balloon-Windkessel hemodynamics.

        Returns a dict with the BOLD time courses for every parameter set,
        a validity mask (``False`` for parameter sets that exploded), and
        optionally the time-averaged SE / SI for E/I-ratio computation.
        """
        print("Simulating MFM 2014 model with simulation parameters:")
        print(
            f"simulate_time: {simulate_time} s, burn_in_time: {burn_in_time} s, "
            f"TR: {TR} s, warm_up_t: {warm_up_t} frames, dt: {self.dt} s"
        )
        N = self.N
        M = self.M
        dt = self.dt

        t_start = 0
        t_end = t_start + burn_in_time + simulate_time
        t_p = torch.arange(t_start, t_end + dt, dt)
        t_len = len(t_p)
        t_inter = int(round(TR / dt))

        z = torch.zeros(N, M)
        f = torch.ones(N, M)
        v_volume = torch.ones(N, M)
        q = torch.ones(N, M)
        bold = torch.zeros(N, M, t_len // t_inter + 1)
        count_bold = 0
        S_E = torch.ones(N, M) * self.S_E_ave
        S_I = torch.ones(N, M) * 0.1433408985

        bold_time_start = time.time()
        print(datetime.datetime.now(), ": Start BOLD calculating...")

        if use_tqdm:
            warm_loop = tqdm(range(warm_up_t), position=0, leave=True)
            main_loop = tqdm(range(t_len), position=0, leave=True)
        else:
            warm_loop = range(warm_up_t)
            main_loop = range(t_len)

        # Warm-up — discard transient dynamics.
        for _ in warm_loop:
            dSE_dt, dSI_dt, _ = self._dynamic_equations(S_E, S_I)
            S_E = S_E + dSE_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            S_I = S_I + dSI_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
        if use_tqdm:
            warm_loop.close()
        print(datetime.datetime.now(), ": End warming up and start main body...")

        S_E_ave = torch.zeros(N, M)
        S_I_ave = torch.zeros(N, M)
        r_E_ave = torch.zeros(N, M)

        for t in main_loop:
            dSE_dt, dSI_dt, r_E = self._dynamic_equations(S_E, S_I)
            S_E = S_E + dSE_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            S_I = S_I + dSI_dt * dt + self.sigma * torch.randn(N, M) * math.sqrt(dt)
            dz_dt, df_dt, dv_dt, dq_dt = self._hemodynamic_equations(
                S_E, z, f, v_volume, q)
            z = z + dz_dt * dt
            f = f + df_dt * dt
            v_volume = v_volume + dv_dt * dt
            q = q + dq_dt * dt

            S_E_ave = S_E_ave + S_E
            S_I_ave = S_I_ave + S_I
            r_E_ave = r_E_ave + r_E

            if (t + 2) % t_inter == 0:
                # Convert hemodynamic state to BOLD via the standard Balloon-Windkessel link.
                bold[:, :, count_bold] = 100 / self.rho * self.V0 * (
                    self.k1 * (1 - q)
                    + self.k2 * (1 - q / v_volume)
                    + self.k3 * (1 - v_volume)
                )
                count_bold += 1
        bold[:, :, count_bold] = 100 / self.rho * self.V0 * (
            self.k1 * (1 - q)
            + self.k2 * (1 - q / v_volume)
            + self.k3 * (1 - v_volume)
        )

        burn_in_bold = int(burn_in_time / TR)
        S_E_ave = S_E_ave / t_len
        S_I_ave = S_I_ave / t_len
        r_E_ave = r_E_ave / t_len

        # Mark exploded / NaN'd parameter sets as invalid (the r_E firing-rate
        # gate was removed in the final paper version — see Sec. 4.2).
        valid_M_mask = torch.ones(M, dtype=torch.bool)
        for i in range(self.M):
            if torch.isnan(r_E_ave[:, i]).any() or torch.isnan(bold[:, i, :]).any():
                valid_M_mask[i] = False

        bold_elapsed = time.time() - bold_time_start
        print("Time using for calculating BOLD signals cost: ", bold_elapsed)
        print("BOLD shape with burn-in: ", bold.shape)
        simulated_bold = bold[:, :, burn_in_bold + 1:]
        print("BOLD shape without burn-in: ", simulated_bold.shape)

        sim_res = {
            "bold": simulated_bold,
            "valid_M_mask": valid_M_mask,
            "r_E_ave": r_E_ave,
        }
        if kwargs.get("need_EI", True):
            sim_res["S_E_ave"] = S_E_ave
            sim_res["S_I_ave"] = S_I_ave
        return sim_res

    # ---------------------------------------------------------------- #
    # Post-processing
    # ---------------------------------------------------------------- #
    def derive_more_dup_sim_res(self, dup_sim_res, **kwargs):
        """
        Aggregate BOLD across the ``param_dup`` duplicates of every parameter set.

        For each parameter set we compute one mean FC matrix, one mean FCD
        histogram, and the time-averaged SE/SI (-> E/I ratio).
        """
        valid_M_mask = dup_sim_res["valid_M_mask"]
        valid_param_indices = get_valid_params_after_dup(dup_sim_res["valid_M_mask"])
        bold_signal = dup_sim_res["bold"]
        r_E_ave = dup_sim_res.get("r_E_ave", None)
        S_E_ave = dup_sim_res.get("S_E_ave", None)
        S_I_ave = dup_sim_res.get("S_I_ave", None)
        save_r_E = r_E_ave is not None
        save_EI = (S_E_ave is not None) and (S_I_ave is not None)

        get_FCD_matrix = kwargs.get("get_FCD_matrix", False)
        get_bold = kwargs.get("get_bold", False)

        num_valid = len(valid_param_indices)
        t_len = bold_signal.shape[3]
        fc_sim = torch.zeros(num_valid, self.N, self.N)
        fcd_hist_sim = torch.zeros(self.fcd_hist_bins, num_valid)
        mean_bold = torch.zeros(self.N, num_valid, t_len)

        if get_FCD_matrix:
            window_num = t_len - self.window_size + 1
            fcd_matrices = torch.zeros(num_valid, window_num, window_num)

        r_E_for_valid_params = torch.zeros(self.N, num_valid) if save_r_E else None
        if save_EI:
            S_E_for_valid_params = torch.zeros(self.N, num_valid)
            S_I_for_valid_params = torch.zeros(self.N, num_valid)
            EI_for_valid_params = torch.zeros(self.N, num_valid)
        else:
            S_E_for_valid_params = S_I_for_valid_params = EI_for_valid_params = None

        for i, idx in enumerate(valid_param_indices):
            mask_this_param = valid_M_mask[idx]
            bold_this_param = bold_signal[:, idx, mask_this_param, :]
            if get_bold:
                mean_bold[:, i] = torch.mean(bold_this_param, dim=1)

            fc_this_param = calc_FC(bold_this_param)
            fc_this_param = torch.mean(fc_this_param, dim=0)
            fcd_mat_this_param, fcd_hist_this_param = calc_FCD(
                bold_this_param, self.window_size, get_FCD_matrix=get_FCD_matrix,
            )
            fcd_hist_this_param = torch.mean(fcd_hist_this_param, dim=1)

            fc_sim[i] = fc_this_param
            fcd_hist_sim[:, i] = fcd_hist_this_param
            if get_FCD_matrix:
                fcd_matrices[i] = fcd_mat_this_param[0]
            if save_r_E:
                r_E_for_valid_params[:, i] = torch.mean(
                    r_E_ave[:, idx, mask_this_param], dim=1)
            if save_EI:
                S_E_ave_this = S_E_ave[:, idx, mask_this_param]
                S_I_ave_this = S_I_ave[:, idx, mask_this_param]
                EI_ave_this = S_E_ave_this / S_I_ave_this
                S_E_for_valid_params[:, i] = torch.mean(S_E_ave_this, dim=1)
                S_I_for_valid_params[:, i] = torch.mean(S_I_ave_this, dim=1)
                EI_for_valid_params[:, i] = torch.mean(EI_ave_this, dim=1)

        save_sim_res = {"valid_param_indices": torch.as_tensor(valid_param_indices)}
        if save_r_E:
            save_sim_res["r_E_for_valid_params"] = r_E_for_valid_params
        if save_EI:
            save_sim_res["S_E_for_valid_params"] = S_E_for_valid_params
            save_sim_res["S_I_for_valid_params"] = S_I_for_valid_params
            save_sim_res["EI_for_valid_params"] = EI_for_valid_params
        if get_FCD_matrix:
            save_sim_res["fc_sim"] = fc_sim
            save_sim_res["fcd_sim"] = fcd_matrices
        if get_bold:
            save_sim_res["bold_TC_sim"] = mean_bold

        other_sim_res = {
            "fc_sim": fc_sim,
            "fcd_hist_sim": fcd_hist_sim,
        }
        return save_sim_res, other_sim_res

    @classmethod
    def calc_losses(cls, save_sim_res, other_sim_res, emp_stats, apply_sort=True):
        """Compute FC corr/L1 + FCD KS losses for each valid parameter set."""
        fc_sim = other_sim_res["fc_sim"]
        fcd_hist_sim = other_sim_res["fcd_hist_sim"]
        emp_fc = emp_stats["emp_fc"]
        emp_fcd_cdf = emp_stats["emp_fcd_cdf"]
        losses = calc_all_loss_from_fc_fcd(fc_sim, fcd_hist_sim, emp_fc, emp_fcd_cdf)

        total_loss = losses["corr_loss"] + losses["l1_loss"] + losses["ks_loss"]
        index_sorted_in_valid = None
        if apply_sort:
            total_loss, index_sorted_in_valid = torch.sort(total_loss, descending=False)
            for key in losses.keys():
                losses[key] = losses[key][index_sorted_in_valid]
            save_sim_res = cls._sort_save_sim_res(save_sim_res, index_sorted_in_valid)
        losses["total_loss"] = total_loss
        save_sim_res.update(losses)
        return save_sim_res

    @classmethod
    def _sort_save_sim_res(cls, save_sim_res, idx):
        """Apply ``idx`` to every dict value (1D along axis 0, 2D along axis 1)."""
        if "valid_param_indices" in save_sim_res:
            save_sim_res["valid_param_indices"] = save_sim_res["valid_param_indices"][idx]
        for k in ("r_E_for_valid_params", "S_E_for_valid_params",
                  "S_I_for_valid_params", "EI_for_valid_params"):
            if k in save_sim_res:
                save_sim_res[k] = save_sim_res[k][:, idx]
        if "fc_sim" in save_sim_res:
            save_sim_res["fc_sim"] = save_sim_res["fc_sim"][idx]
        if "fcd_sim" in save_sim_res:
            save_sim_res["fcd_sim"] = save_sim_res["fcd_sim"][idx]
        if "bold_TC_sim" in save_sim_res:
            save_sim_res["bold_TC_sim"] = save_sim_res["bold_TC_sim"][:, idx]
        return save_sim_res

    # ---------------------------------------------------------------- #
    # FIC equations
    # ---------------------------------------------------------------- #
    def _dynamic_equations(self, S_E_t, S_I_t):
        """Deco 2014 FIC equations — returns dSE/dt, dSI/dt and the firing rate r_E."""
        I_E_t = (
            self.W_E * self.I_0
            + self.w_EE * self.J_NMDA * S_E_t
            + self.G * self.J_NMDA * torch.matmul(self.sc_euler, S_E_t)
            - self.w_IE * S_I_t
        )
        I_I_t = self.W_I * self.I_0 + self.w_EI * self.J_NMDA * S_E_t - S_I_t
        r_E = (self.a_E * I_E_t - self.b_E) / (
            1 - torch.exp(-self.d_E * (self.a_E * I_E_t - self.b_E))
        )
        r_I = (self.a_I * I_I_t - self.b_I) / (
            1 - torch.exp(-self.d_I * (self.a_I * I_I_t - self.b_I))
        )
        dSE_dt = -S_E_t / self.tau_E + (1 - S_E_t) * self.gamma_kin * r_E
        dSI_dt = -S_I_t / self.tau_I + r_I
        return dSE_dt, dSI_dt, r_E

    def _hemodynamic_equations(self, S_t, z_t, f_t, v_t, q_t):
        """Balloon-Windkessel hemodynamic forward model."""
        dz_dt = S_t - self.kappa * z_t - self.gamma_hemo * (f_t - 1)
        df_dt = z_t
        dv_dt = (f_t - v_t ** (1 / self.alpha)) / self.tau
        dq_dt = (
            f_t / self.rho * (1 - (1 - self.rho) ** (1 / f_t))
            - q_t * v_t ** (1 / self.alpha - 1)
        ) / self.tau
        return dz_dt, df_dt, dv_dt, dq_dt

    def _solve_I_I(self, I_I_ave, S_E_ave, w_EI_m):
        if torch.cuda.is_available():
            w_EI_m = w_EI_m.cpu().numpy()
        phi_I_I_ave = (self.a_I * I_I_ave - self.b_I) / (
            1 - np.exp(-self.d_I * (self.a_I * I_I_ave - self.b_I))
        )
        return (self.W_I * self.I_0
                + self.J_NMDA * w_EI_m * S_E_ave
                - phi_I_I_ave * self.tau_I
                - I_I_ave)

    def _solve_S_E_ave(self, S_E_ave):
        return S_E_ave / (self.tau_E * self.gamma_kin * (1 - S_E_ave)) - self.r_E

    def _solve_I_E_ave(self, I_E_ave):
        tmp = self.a_E * I_E_ave - self.b_E
        return tmp / (1 - np.exp(-self.d_E * tmp)) - self.r_E


# Registry consulted by :func:`DELSSOME_indiv.optimizers.CBIG_cmaes.ModelOptimizer` to
# locate the dynamic-model class by its short name.
AVAIL_MODELS = {"pfic": pFIC}
