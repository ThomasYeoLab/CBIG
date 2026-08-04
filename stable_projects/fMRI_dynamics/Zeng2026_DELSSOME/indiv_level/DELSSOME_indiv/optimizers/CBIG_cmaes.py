# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
CMA-ES parameter inversion of the FIC model — Euler or DELSSOME.

Highlights:

* :class:`ModelOptimizer` is the abstract base: it owns the dynamic model,
  the empirical statistics for one subject (or group), and the dt for the
  current phase (train / val / test — different dt, see ``config.ini``).

* :class:`CmaesTrainer` implements the original CMA-ES of Hansen (2006). On
  each epoch ``k`` it:
    1. Samples ``param_sets=100`` candidate parameter-coefficient vectors from
       :math:`\\mathcal{N}(m_k, \\sigma_k^2 C_k)`, rejecting out-of-range or
       low-variance ``wEI`` profiles (pFIC.is_parameter_valid).
    2. Evaluates each candidate using EITHER Euler integration (epoch in
       ``euler_epoch_range``) OR DELSSOME (epoch in ``dl_epoch_range``).
    3. Picks the best ``select_param_sets=10`` candidates, dumps the per-epoch
       parameter dict (used downstream to train DELSSOME), and uses them to
       update ``m_{k+1}``, ``\\sigma_{k+1}``, ``C_{k+1}``.

* :class:`ModelValidator` and :class:`ModelTester` re-simulate the best
  candidates from the previous phase with a smaller ``dt`` (more accurate).

* :func:`train_helper` retries the whole optimisation a few times if a phase
  collapses (e.g. CMA-ES could not find enough valid candidates).
"""

import time
import shutil
import os
import datetime

import numpy as np
import pandas as pd
import torch
import torch.distributions as td

from DELSSOME_indiv.models.CBIG_dynamic_model import (
    DynamicModel, reshape_dup_sim_res, AVAIL_MODELS,
)
from DELSSOME_indiv.CBIG_constants import PREV_PHASE
from DELSSOME_indiv.CBIG_file_utils import (
    get_best_params_file_path, get_train_epoch_file_path,
)
from DELSSOME_indiv.CBIG_misc_utils import (
    parameterize_myelin_rsfc,  # noqa: F401
    set_torch_default,
)


def csv2tensor(csv_path):
    """Load a 2D CSV as a torch double tensor (no header / index)."""
    content = pd.read_csv(csv_path, header=None, index_col=False)
    return torch.as_tensor(np.array(content), dtype=torch.double)


def select_best_parameter_from_savedict(save_dict):
    """Return the lowest-total-loss parameter from a CMA-ES save dict."""
    if "corr_loss" in save_dict:
        total_loss = save_dict["corr_loss"] + save_dict["L1_loss"] + save_dict["ks_loss"]
    elif "pred_loss" in save_dict:
        total_loss = torch.sum(save_dict["pred_loss"], dim=1)
    elif "FC_FCD_loss" in save_dict:
        total_loss = torch.sum(save_dict["FC_FCD_loss"], dim=1)
    else:
        raise Exception("Save dict has neither 'corr_loss' nor 'pred_loss'.")

    index_in_valid = torch.argmin(total_loss)
    index_min = save_dict["valid_param_indices"][index_in_valid]
    if "parameter" in save_dict:
        return save_dict["parameter"][:, index_min]
    if "param_coef" in save_dict:
        return parameterize_myelin_rsfc(save_dict["param_coef"][:, index_min])
    raise Exception("Save dict has neither 'parameter' nor 'param_coef'.")


# --------------------------------------------------------------------------- #
# Retry wrapper
# --------------------------------------------------------------------------- #
def train_helper(config, emp_stats, train_save_dir, num_epochs,
                 dl_epoch_range, euler_epoch_range, opportunities,
                 next_epoch, seed=None, other_parameterization=None):
    """Try CMA-ES up to ``opportunities`` times before giving up.

    Returns the trainer's exit state:
      * 0 — success;
      * 1 — could not initialise CMA-ES (too few valid candidates at startup);
      * 2 — CMA-ES broke in the middle (re-run from scratch).
    """
    state = 1
    for i in range(opportunities):
        mfm_trainer = CmaesTrainer(
            config=config, emp_stats=emp_stats,
            train_save_dir=train_save_dir, num_epochs=num_epochs,
            dl_epoch_range=dl_epoch_range, euler_epoch_range=euler_epoch_range,
            other_parameterization=other_parameterization,
        )
        state = mfm_trainer.train(seed=seed, next_epoch=next_epoch)
        if state == 0:
            break
        if state == 1:
            if i == opportunities - 1:
                print(f"Having tried for {opportunities} times, but still not enough "
                      "parameter to initialise CMA-ES.")
            print("Failed to initialise CMA-ES. Restarting ...")
        elif state == 2:
            print("CMA-ES broke during middle epochs. Restarting...")
    print(f"In total {i + 1} random seeds.")
    print("The End.")
    return state


# --------------------------------------------------------------------------- #
# Base class
# --------------------------------------------------------------------------- #
class ModelOptimizer:
    """
    Common state for Trainer / Validator / Tester.

    Builds the dynamic-model class from ``config['Dynamic Model Constants']
    ['dynamic_model_name']`` (``pfic`` here), reads simulation parameters,
    and stores empirical SC/FC/FCD that the per-epoch logic will need.
    """

    def __init__(self, config, phase, emp_stats, prev_phase_best_params_path,
                 curr_phase_save_dir):
        assert phase in ["train", "val", "test"], (
            "phase must be one of ['train', 'val', 'test']")
        self.device, self.dtype = set_torch_default()
        self.config = config
        self.phase = phase
        self.curr_phase_save_dir = curr_phase_save_dir
        self.prev_phase_best_params_path = prev_phase_best_params_path
        os.makedirs(self.curr_phase_save_dir, exist_ok=True)

        self.dynamic_model_name = config["Dynamic Model Constants"]["dynamic_model_name"]
        self.dynamic_model_cls: DynamicModel = AVAIL_MODELS[self.dynamic_model_name]

        # Dataset / simulation parameters
        dp = config["Dataset Parameters"]
        self.simulate_time = float(dp["simulate_time"])
        self.burn_in_time = float(dp["burn_in_time"])
        self.TR = float(dp["TR"])
        self.window_size = int(dp["window_size"])

        sp = config["Simulating Parameters"]
        self.N = int(sp["n_ROI"])
        self.parameter_dim = self.dynamic_model_cls.get_parameter_dim(self.N)
        self.param_dup = int(sp["param_dup"])
        self.warm_up_t = int(sp["warm_up_t"])
        self.dt = float(sp[f"dt_{phase}"])

        # Empirical statistics
        self.emp_stats = emp_stats
        self.sc_euler = emp_stats["sc_euler"]
        self.fc_euler = emp_stats["emp_fc"]
        self.fcd_cdf_euler = emp_stats["emp_fcd_cdf"]

        if self.phase == "train":
            self.sc_dl = emp_stats["sc_dl"]
            self.fc_dl = emp_stats["fc_dl"]
            self.fcd_dl = emp_stats["fcd_dl"]
        print(self.post_init_message)

    @property
    def post_init_message(self):
        return (f"Successfully init ModelOptimizer at phase {self.phase}. "
                f"The results will be saved in {self.curr_phase_save_dir}")

    @property
    def prev_phase(self):
        return PREV_PHASE[self.phase]

    # ----- Euler-based candidate evaluation -------------------------------- #
    def sim_param_with_dup(self, param_vectors, reshape_res=True):
        """Run Euler integration for ``param_vectors`` with each parameter set
        duplicated ``param_dup`` times (different noise realisations)."""
        input_param_sets = param_vectors.shape[1]
        parameter_repeat = param_vectors.repeat(1, self.param_dup)
        dynamic_model = self.dynamic_model_cls(
            self.config, parameter_repeat, self.sc_euler, dt=self.dt
        )
        sim_res = dynamic_model.simulate(
            simulate_time=self.simulate_time,
            burn_in_time=self.burn_in_time,
            TR=self.TR,
            warm_up_t=self.warm_up_t,
        )
        if reshape_res:
            sim_res = reshape_dup_sim_res(sim_res, input_param_sets, self.param_dup)
        return dynamic_model, sim_res

    def sim_n_get_losses(self, param_vectors, sort_param_vectors=True,
                         get_FCD_matrix=False, get_bold=False):
        """Run Euler integration and compute FC+FCD losses for every param set."""
        dynamic_model, dup_sim_res = self.sim_param_with_dup(
            param_vectors, reshape_res=True)
        save_sim_res, other_sim_res = dynamic_model.derive_more_dup_sim_res(
            dup_sim_res, get_FCD_matrix=get_FCD_matrix, get_bold=get_bold)
        save_sim_res = dynamic_model.calc_losses(
            save_sim_res, other_sim_res, self.emp_stats, apply_sort=True)
        if len(save_sim_res["valid_param_indices"]) == 0:
            raise Exception(
                "No valid parameter vectors found (likely due to all parameters "
                "being out-of-range)!"
            )
        save_sim_res["parameter"] = (
            param_vectors[:, save_sim_res["valid_param_indices"]]
            if sort_param_vectors else param_vectors
        )
        save_sim_res["are_params_sorted"] = sort_param_vectors
        return save_sim_res


# --------------------------------------------------------------------------- #
# CMA-ES Trainer
# --------------------------------------------------------------------------- #
class CmaesTrainer(ModelOptimizer):
    """
    Optimise FIC parameters with CMA-ES (Hansen, 2006).

    On each epoch ``k``:
      * Sample 100 candidate parameter-coefficient vectors from
        ``Normal(m_k, sigma_k^2 * C_k)``, retrying until 100 of them pass
        :meth:`pFIC.is_parameter_valid`.
      * Evaluate them with Euler integration OR DELSSOME (depending on whether
        ``k`` falls in ``euler_epoch_range`` or ``dl_epoch_range``).
      * Pick the top-10 by total loss, save a per-epoch dump, and update
        the CMA-ES distribution from those 10 winners.
    """

    def __init__(self, config, emp_stats, train_save_dir,
                 num_epochs=100, dl_epoch_range=(),
                 euler_epoch_range=tuple(range(100)),
                 other_parameterization=None):
        super().__init__(config, "train", emp_stats, train_save_dir, train_save_dir)

        # CMA-ES hyper-parameters
        self.param_sets = 100             # 100 candidates per epoch
        self.select_param_sets = 10       # top-10 update the distribution
        self.min_select_param_sets = 2    # absolute minimum to not collapse

        # Parameterisation
        self.myelin = emp_stats["myelin"]
        self.rsfc_gradient = emp_stats["rsfc_gradient"]
        self.other_parameterization = other_parameterization
        (self.param_coef_dim,
         self.concat_mat,
         self.pinv_concat_mat) = self.dynamic_model_cls.get_parameterization(
            emp_stats=self.emp_stats,
            other_parameterization=other_parameterization,
        )

        # Per-epoch routing: which epochs use Euler, which use DELSSOME?
        self.num_epochs = int(num_epochs)
        self.dl_epoch_range = np.array(dl_epoch_range, dtype=int)
        self.euler_epoch_range = np.array(euler_epoch_range, dtype=int)
        self.epoch_range = np.concatenate((self.dl_epoch_range, self.euler_epoch_range))

        if len(self.dl_epoch_range) > 0:
            self.dl_models = self.dynamic_model_cls.load_dl_models(
                config, self.emp_stats, self.device
            )

    # ----- per-epoch parameter dump --------------------------------------- #
    def select_valid_params(self, save_dict, epoch, is_saving=True):
        valid_param_indices = save_dict["valid_param_indices"]
        total_num_valid = len(valid_param_indices)
        if total_num_valid < self.min_select_param_sets:
            print(f"Too few valid param vectors to continue "
                  f"(only {total_num_valid} parameters)!")
            return None, None
        if total_num_valid < self.select_param_sets:
            used_num_valid = self.min_select_param_sets
            print(f"Although we do not have enough valid param vectors "
                  f"({total_num_valid}), we will still continue using "
                  f"{used_num_valid} param vectors.")
        else:
            used_num_valid = self.select_param_sets

        if is_saving:
            print("Start saving parameter, valid_param_indices and losses...")
            train_file_path = get_train_epoch_file_path(self.curr_phase_save_dir, epoch)
            torch.save(save_dict, train_file_path)
            print("Successfully saved params and losses to:", train_file_path)

        return (
            save_dict["total_loss"][:used_num_valid],
            valid_param_indices[:used_num_valid],
        )

    # ----- main optimisation loop ----------------------------------------- #
    def train(self, seed=None, next_epoch=0):
        if next_epoch >= self.num_epochs:
            raise Exception("You do not need any more epoch.")

        N = self.N
        self.search_range = self.dynamic_model_cls.get_search_range(N, self.emp_stats)

        if next_epoch == 0:
            # Start from scratch: clear the save dir.
            if os.path.exists(self.curr_phase_save_dir):
                print(f"Deleting {self.curr_phase_save_dir} ...")
                shutil.rmtree(self.curr_phase_save_dir)
            os.makedirs(self.curr_phase_save_dir, exist_ok=True)

            if seed is None:
                seed = np.random.randint(0, 1000000000000)
            rand_ge = torch.manual_seed(seed)

            (m_k, sigma_k, cov_k, p_sigma_k, p_c_k) = self.dynamic_model_cls.init_CMA_ES(
                self.search_range, self.pinv_concat_mat, self.param_coef_dim
            )
            print("Training parameters have been initialised...")

        elif next_epoch in self.euler_epoch_range:
            previous_final_state_path = os.path.join(
                self.curr_phase_save_dir, "final_state.pth")
            if not os.path.exists(previous_final_state_path):
                raise Exception("Previous final state path doesn't exist.")
            final_dict = torch.load(previous_final_state_path, map_location=self.device)
            seed = final_dict["seed"]
            random_state = final_dict["random_state"].to(
                torch.ByteTensor(device=torch.device("cpu")))
            m_k = final_dict["m"]
            sigma_k = final_dict["sigma"]
            cov_k = final_dict["cov"]
            p_sigma_k = final_dict["p_sigma"]
            p_c_k = final_dict["p_c"]

            rand_ge = torch.manual_seed(seed)
            rand_ge.set_state(random_state)
            print("Successfully loaded previous parameters. Will start next epochs.")
        else:
            raise Exception("Argument next_epoch invalid.")

        print(" -- Start CMA-ES optimisation --")
        for k in range(next_epoch, self.num_epochs):
            print(f"Epoch: [{k + 1}/{self.num_epochs}]")
            if k in self.epoch_range:
                start_time = time.time()
                epoch_res = self._train_one_epoch(
                    k, m_k, sigma_k, cov_k, p_sigma_k, p_c_k)
                print(f"Epoch {k} elapsed time: {time.time() - start_time}")
                if epoch_res is None:
                    # State 1 if the very first epoch failed (init problem),
                    # state 2 otherwise (mid-optimisation collapse).
                    return (
                        1 if k in self.dl_epoch_range
                        or k == self.euler_epoch_range[0]
                        else 2
                    )
                m_k, sigma_k, cov_k, p_sigma_k, p_c_k = epoch_res

        print(" -- CMA-ES completed --")
        final_dict = {
            "seed": seed,
            "random_state": rand_ge.get_state(),
            "m": m_k, "sigma": sigma_k, "cov": cov_k,
            "p_sigma": p_sigma_k, "p_c": p_c_k,
            "epoch": k,
        }
        torch.save(final_dict,
                   os.path.join(self.curr_phase_save_dir, "final_state.pth"))
        print("Random seed, random state and final CMA-ES parameters saved.")
        return 0

    def _train_one_epoch(self, k, m_k, sigma_k, cov_k, p_sigma_k, p_c_k):
        cov_for_gaussian = sigma_k ** 2 * cov_k
        if torch.isnan(cov_for_gaussian).any():
            print("cov_for_gaussian contains nan!")
            print("sigma_k: ", sigma_k)
            print("cov_k: ", cov_k)
            return None
        param_coef_k, parameter_k = self._sample_valid_parameters(
            m_k, cov_for_gaussian, self.search_range)
        if param_coef_k is None or parameter_k is None:
            print("Sampling failed!")
            return None
        print(f"Successfully sampled {param_coef_k.shape[1]} parameters.", flush=True)

        if k in self.dl_epoch_range:
            # DELSSOME evaluation (Section 2.3 of the paper).
            save_sim_res = self.dynamic_model_cls.sim_n_get_losses_dl(
                dl_models=self.dl_models, param_vectors=parameter_k,
                emp_stats=self.emp_stats, sort_param_vectors=False,
            )
        elif k in self.euler_epoch_range:
            save_sim_res = self.sim_n_get_losses(parameter_k, sort_param_vectors=False)
        else:
            raise Exception("Invalid epoch. Break.")

        loss_k, index_k = self.select_valid_params(save_sim_res, k)
        if loss_k is None:
            print(f"No loss returned, iteration {k} ends.")
            return None
        select_params = param_coef_k[:, index_k]

        return self._update_CMA_ES(
            select_params, loss_k, m_k, sigma_k, cov_k, p_sigma_k, p_c_k, k
        )

    def _update_CMA_ES(self, select_params, loss_k, m_k, sigma_k, cov_k,
                       p_sigma_k, p_c_k, k):
        """Standard CMA-ES update — see https://en.wikipedia.org/wiki/CMA-ES."""
        select_param_sets = select_params.shape[1]
        loss_inverse = 1 / loss_k
        weights = loss_inverse / torch.sum(loss_inverse)
        m_kp1 = torch.matmul(select_params, weights)
        mueff = 1 / torch.sum(weights ** 2)

        # Eigendecomposition for C^{-1/2}.
        Lambda, V = torch.linalg.eigh(cov_k)
        Lambda = torch.sqrt(Lambda)
        inv_sqrt_cov = torch.matmul(V, torch.matmul(torch.diag(Lambda ** -1), V.T))

        # Standard CMA-ES coefficients.
        c_sigma = (mueff + 2) / (self.param_coef_dim + mueff + 5)
        c_c = (4 + mueff / self.param_coef_dim) / (
            self.param_coef_dim + 4 + 2 * mueff / self.param_coef_dim
        )
        c_1 = 2 / ((self.param_coef_dim + 1.3) ** 2 + mueff)
        c_mu = min(1 - c_1,
                   2 * (mueff - 2 + 1 / mueff) /
                   ((self.param_coef_dim + 2) ** 2 + mueff))
        d_sigma = 1 + 2 * max(
            0, torch.sqrt((mueff - 1) / (self.param_coef_dim + 1)) - 1
        ) + c_sigma
        expected_value = self.param_coef_dim ** 0.5 * (
            1 - 1 / (4 * self.param_coef_dim)
            + 1 / (21 * self.param_coef_dim ** 2)
        )

        # Update evolution paths.
        p_sigma_kp1 = (1 - c_sigma) * p_sigma_k + torch.sqrt(
            c_sigma * (2 - c_sigma) * mueff
        ) * torch.matmul(inv_sqrt_cov, (m_kp1 - m_k).unsqueeze(1) / sigma_k)
        indicator = (
            torch.linalg.norm(p_sigma_kp1)
            / torch.sqrt(1 - (1 - c_sigma) ** (2 * k))
            / expected_value
            < (1.4 + 2 / (self.param_coef_dim + 1))
        ) * 1
        p_c_kp1 = (1 - c_c) * p_c_k + indicator * torch.sqrt(
            c_c * (2 - c_c) * mueff
        ) * (m_kp1 - m_k).unsqueeze(1) / sigma_k

        # Adapt C and sigma.
        artmp = (1 / sigma_k) * (
            select_params - torch.tile(m_k, [select_param_sets, 1]).T
        )
        cov_kp1 = (
            (1 - c_1 - c_mu) * cov_k
            + c_1 * (torch.matmul(p_c_kp1, p_c_kp1.T)
                     + (1 - indicator) * c_c * (2 - c_c) * cov_k)
            + c_mu * torch.matmul(artmp,
                                  torch.matmul(torch.diag(weights), artmp.T))
        )
        sigma_kp1 = sigma_k * torch.exp(
            (c_sigma / d_sigma)
            * (torch.linalg.norm(p_sigma_kp1) / expected_value - 1)
        )
        return m_kp1, sigma_kp1, cov_kp1, p_sigma_kp1, p_c_kp1

    def _sample_valid_parameters(self, mean, cov, search_range):
        """Sample 100 parameter vectors that pass ``is_parameter_valid``."""
        multivariate_normal = td.MultivariateNormal(mean, cov)
        valid_count = 0
        total_count = 0
        total_threshold = 20000 * self.param_sets
        sampled_params = torch.zeros(self.param_coef_dim, self.param_sets)
        sampled_parameters = torch.zeros(self.parameter_dim, self.param_sets)

        while valid_count < self.param_sets:
            sampled_params[:, valid_count] = multivariate_normal.sample()
            sampled_parameters[:, valid_count] = self.dynamic_model_cls.get_parameters(
                sampled_params[:, valid_count], self.concat_mat
            )
            if self.dynamic_model_cls.is_parameter_valid(
                self.config, self.emp_stats,
                sampled_parameters[:, valid_count], search_range,
            ):
                valid_count += 1
            total_count += 1
            if total_count >= total_threshold:
                print(f"Not enough valid sampled parameters! "
                      f"Only sample {valid_count} parameters!")
                return None, None
        return sampled_params, sampled_parameters


# --------------------------------------------------------------------------- #
# Validation / Test
# --------------------------------------------------------------------------- #
class ModelValidator(ModelOptimizer):
    """Re-simulate ``best_from_<prev_phase>`` candidates with the validation dt."""

    def __init__(self, config, emp_stats, prev_phase_best_params_path,
                 curr_phase_save_dir, phase="val"):
        super().__init__(config, phase, emp_stats,
                         prev_phase_best_params_path, curr_phase_save_dir)

    def validate(self, use_top_k=None, seed=None):
        print(f" -- Start {self.phase} phase -- ")
        if seed is None:
            seed = np.random.randint(0, 1000000000000)
        torch.manual_seed(seed)

        best_from_prev_phase = torch.load(
            self.prev_phase_best_params_path, map_location=self.device
        )
        param_vectors = best_from_prev_phase["parameter"]
        if use_top_k is not None:
            param_vectors = param_vectors[:, :use_top_k]
        save_dict = self.sim_n_get_losses(param_vectors)
        save_dict["seed"] = seed

        # Carry over the previous-phase losses for traceability.
        valid_param_indices = save_dict["valid_param_indices"]
        for key, value in best_from_prev_phase.items():
            if key.endswith("loss"):
                value = value[valid_param_indices]
                if key.startswith("train"):
                    save_dict[key] = value
                else:
                    save_dict[f"{self.prev_phase}_{key}"] = value

        print(datetime.datetime.now(), "Start saving results...")
        sim_on_curr_phase_file_path = get_best_params_file_path(
            self.phase, self.curr_phase_save_dir
        )
        torch.save(save_dict, sim_on_curr_phase_file_path)
        print("Successfully saved to:", sim_on_curr_phase_file_path)
        print(datetime.datetime.now(), f" -- Done {self.phase} phase -- ")
        return 0


class ModelTester(ModelValidator):
    """Test = same logic as validation but uses only the single best parameter."""

    def __init__(self, config, emp_stats, prev_phase_best_params_path,
                 curr_phase_save_dir):
        super().__init__(
            config=config, emp_stats=emp_stats,
            prev_phase_best_params_path=prev_phase_best_params_path,
            curr_phase_save_dir=curr_phase_save_dir, phase="test",
        )
        self.test_param_sets = 1

    def test(self, use_top_k=None, seed=None):
        if use_top_k is None:
            use_top_k = self.test_param_sets
        return self.validate(use_top_k=use_top_k, seed=seed)
