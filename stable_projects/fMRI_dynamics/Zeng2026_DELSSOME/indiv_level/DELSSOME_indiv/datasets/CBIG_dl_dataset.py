# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
PyTorch / Lightning DataModules for training the DELSSOME cost predictor.

The data flow is:

* For each dataset listed in ``ds_names`` and each phase (train/val/test):
   * Read the subject list from
     ``LOG_DIR/<ds>/train_indiv_DELSSOME-pfic/<phase>/sub_list.txt``.
   * Load the (per-subject) empirical SC + FC + FCD-PDF via
     :func:`DELSSOME_indiv.CBIG_file_utils.get_stats(with_dl_stats=True)`.
   * For every Euler-CMA-ES training epoch ``k``, open
     ``LOG_DIR/<ds>/indiv_Euler-pfic/train/trial1/seed1/<sub_id>/param_save_epoch<k>.pth``
     and emit one training sample per "child" (candidate parameter set):
        - **input**:  ``(ds_idx, params, sc_euler, emp_fc, fcd_dl, G)`` for the
          multimodal predictor or ``(ds_idx, param_vector, sc_dl, fc_dl, fcd_dl)``
          for the MLP-only baseline.
        - **target**: ``[corr_loss, l1_loss, ks_loss]``.

* Training samples from different datasets and epochs are interleaved so each
  mini-batch sees a mix of contexts.

The resulting in-memory tensors are cached under ``TMP_DIR/cache`` (override with
the ``DELSSOME_TMP_DIR`` env var, so a read-only data mount still works); the
cache key embeds the dataset list and ``n_sample_per_phase``.
"""

import os
import warnings
import torch
import numpy as np
import matplotlib.pyplot as plt
from lightning import LightningDataModule
from torch.utils.data import DataLoader, Dataset

from DELSSOME_indiv.CBIG_constants import ALL_TRAIN_DS, LOG_DIR, TMP_DIR, get_ds_idx
from DELSSOME_indiv.CBIG_file_utils import (
    get_run_dir,
    get_stats,
    get_train_epoch_file_path,
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def n_sample_per_phase_to_list(n_sample_per_phase, all_sub_or_group_list, phase):
    """Slice an existing subject list to the number of samples requested."""
    n_sample = 0
    for d in n_sample_per_phase:
        for p, n in d.items():
            if p == phase:
                n_sample += n
    if n_sample > len(all_sub_or_group_list):
        warnings.warn(
            f"# samples per phase ({n_sample}) > # subjects/groups "
            f"({len(all_sub_or_group_list)}). Using all available subjects."
        )
    if n_sample == 0:
        raise ValueError(f"Number of samples per phase ({n_sample}) is 0")
    return all_sub_or_group_list[:n_sample]


def n_sample_per_phase_to_file_name(n_sample_per_phase):
    """Short string describing the (train, val, test) sample sizes."""
    file_name = ""
    for d in n_sample_per_phase:
        for p, n in d.items():
            file_name += f"{p[:2]}{n}"
    return file_name


# --------------------------------------------------------------------------- #
# Base classes
# --------------------------------------------------------------------------- #
class ParamDataset(Dataset):
    """Concatenate per-(subject, epoch) param dumps into one flat dataset.

    Caching: the materialised ``(data, targets)`` is dumped to disk on first use
    and reloaded on every subsequent call.  Cache key includes the dataset list
    and the ``n_sample_per_phase`` spec.
    """

    def __init__(self,
                 phase,
                 ds_names=ALL_TRAIN_DS,
                 group_type=None,
                 n_epochs=100,
                 overwrite_cache=False,
                 default_dtype=torch.float32,
                 used_n_children=1,
                 n_sample_per_phase=[{"train": 1, "val": 1, "test": 1}]):
        self.used_n_children = used_n_children
        self.phase = phase
        self.ds_names = ds_names
        self.group_type = group_type
        self.default_dtype = default_dtype
        self.n_sample_per_phase = n_sample_per_phase
        self.data = []
        self.targets = []

        if not overwrite_cache and os.path.exists(self.cache_path):
            self.load_from_cache()
        else:
            self._load_and_process_data(n_epochs)
            self.save_to_cache()

    # ----- subclass hooks --------------------------------------------------- #
    def _get_list_file_path(self, ds):
        raise NotImplementedError

    def _get_param_folder(self, ds, group_or_sub_id):
        raise NotImplementedError

    def param_dict_to_data(self, ds_idx, param_dict, used_n_children, indiv_stats):
        raise NotImplementedError

    # ----- core loading logic ---------------------------------------------- #
    def _load_and_process_data(self, n_epochs):
        all_data, all_targets = self._collect_all_data(n_epochs)
        self._interleave_data(all_data, all_targets, n_epochs)

    def _collect_all_data(self, n_epochs):
        all_data, all_targets = {}, {}
        for ds in self.ds_names:
            sub_or_group_list = self._get_sub_or_group_list(ds)
            all_data[ds], all_targets[ds] = {}, {}
            for group_or_sub_id in sub_or_group_list:
                stats = get_stats(
                    ds_name=ds, group_or_sub_id=group_or_sub_id,
                    group_type=self.group_type, atlas="DK68",
                    with_myelin_gradient=False, with_dl_stats=True,
                    dtype=self.default_dtype,
                )
                ds_idx = torch.tensor(get_ds_idx(ds))
                for epoch_idx in range(n_epochs):
                    all_data[ds].setdefault(epoch_idx, [])
                    all_targets[ds].setdefault(epoch_idx, [])
                    param_dict = self._load_param_dict(ds, group_or_sub_id, epoch_idx)
                    if param_dict is None:
                        continue
                    data, targets = self.param_dict_to_data(
                        ds_idx, param_dict, self.used_n_children, stats
                    )
                    all_data[ds][epoch_idx].extend(data)
                    all_targets[ds][epoch_idx].extend(targets)
        return all_data, all_targets

    def _get_sub_or_group_list(self, ds):
        list_file_path = self._get_list_file_path(ds)
        with open(list_file_path, "r") as f:
            sub_or_group_list = [sub.strip() for sub in f.readlines()]
        return n_sample_per_phase_to_list(
            self.n_sample_per_phase, sub_or_group_list, self.phase
        )

    def _load_param_dict(self, ds, group_or_sub_id, epoch_idx):
        param_folder = self._get_param_folder(ds, group_or_sub_id)
        param_dict_path = get_train_epoch_file_path(param_folder, epoch_idx)
        if not os.path.exists(param_dict_path):
            print(f"Warning: {param_dict_path} does not exist.")
            return None
        return torch.load(param_dict_path, map_location="cpu")

    def _interleave_data(self, all_data, all_targets, n_epochs):
        # Round-robin over (epoch, dataset) so a mini-batch mixes datasets.
        for epoch_idx in range(n_epochs):
            for ds in self.ds_names:
                if ds in all_data and epoch_idx in all_data[ds]:
                    self.data.extend(all_data[ds][epoch_idx])
                    self.targets.extend(all_targets[ds][epoch_idx])

    # ----- caching --------------------------------------------------------- #
    @property
    def cache_path_prefix(self):
        prefix = f"{self.phase}-"
        if self.ds_names != ALL_TRAIN_DS:
            prefix = f"{'_'.join(self.ds_names)}-{prefix}"
        return prefix

    @property
    def cache_path_suffix(self):
        return (f"-n_child{self.used_n_children}-"
                f"{n_sample_per_phase_to_file_name(self.n_sample_per_phase)}")

    @property
    def cache_path(self):
        return os.path.join(
            TMP_DIR, "cache",
            f"{self.cache_path_prefix}ds{self.cache_path_suffix}.pth",
        )

    def save_to_cache(self):
        os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
        torch.save((self.data, self.targets), self.cache_path)

    def load_from_cache(self):
        self.data, self.targets = torch.load(self.cache_path)
        print(f"Loaded data from {self.cache_path}")

    # ----- Dataset API ----------------------------------------------------- #
    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        return self.data[index], self.targets[index]


class DataModule(LightningDataModule):
    """Generic Lightning DataModule wrapper around a ``ParamDataset`` subclass."""

    def __init__(self, is_trial_run=False, ds_names=ALL_TRAIN_DS,
                 group_type=None,
                 n_sample_per_phase=[{"train": 1, "val": 1, "test": 1}],
                 batch_size=256, num_workers=0):
        super().__init__()
        if is_trial_run:
            print("Running trial run!")
            n_sample_per_phase = [{"train": 1, "val": 1, "test": 1}]  # noqa: F841 (captured by save_hyperparameters)
        self.pin_memory = num_workers > 0
        self.save_hyperparameters()

    @property
    def ds_class(self):
        raise NotImplementedError

    def setup(self, stage):
        used_n_children = 1 if self.hparams.is_trial_run else 100
        kwargs = dict(
            ds_names=self.hparams.ds_names,
            group_type=self.hparams.group_type,
            used_n_children=used_n_children,
            n_sample_per_phase=self.hparams.n_sample_per_phase,
        )
        if stage == "fit":
            self.train_ds = self.ds_class(phase="train", **kwargs)
            self.val_ds = self.ds_class(phase="val", **kwargs)
        if stage == "test":
            self.test_ds = self.ds_class(phase="test", **kwargs)
        if stage == "predict":
            self.predict_ds = self.ds_class(phase="test", **kwargs)

    def train_dataloader(self):
        return DataLoader(self.train_ds, batch_size=self.hparams.batch_size,
                          num_workers=self.hparams.num_workers,
                          pin_memory=self.pin_memory, shuffle=True)

    def val_dataloader(self):
        return DataLoader(self.val_ds, batch_size=self.hparams.batch_size,
                          num_workers=self.hparams.num_workers,
                          pin_memory=self.pin_memory)

    def test_dataloader(self):
        return DataLoader(self.test_ds, batch_size=self.hparams.batch_size,
                          num_workers=self.hparams.num_workers,
                          pin_memory=self.pin_memory)


# --------------------------------------------------------------------------- #
# MLP-only predictor (params, SC vec, FC vec, FCD pdf) -> 3 losses
# --------------------------------------------------------------------------- #
class PredictorDataset(ParamDataset):

    @property
    def cache_path_prefix(self):
        return f"pred-{super().cache_path_prefix}"

    def param_dict_to_data(self, ds_idx, param_dict, used_n_children, emp_stats):
        data, targets = [], []
        sc_dl = emp_stats["sc_dl"]
        fc_dl = emp_stats["fc_dl"]
        fcd_dl = emp_stats["fcd_dl"]
        for i, child_idx in enumerate(
                param_dict["valid_param_indices"][:used_n_children]):
            param_vector = param_dict["parameter"][:, child_idx].to(self.default_dtype)
            data.append((ds_idx, param_vector, sc_dl, fc_dl, fcd_dl))
            corr_loss = param_dict["corr_loss"][i]
            l1_loss = param_dict["l1_loss"][i]
            ks_loss = param_dict["ks_loss"][i]
            targets.append(torch.tensor([corr_loss, l1_loss, ks_loss],
                                        dtype=self.default_dtype))
        return data, targets

    def plot_dist(self):
        targets = np.array(self.targets)
        loss_types = ["corr_loss", "l1_loss", "ks_loss"]
        colors = ["blue", "green", "orange", "purple"]
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.flatten()
        for i, (loss_type, color) in enumerate(zip(loss_types, colors)):
            axes[i].hist(targets[:, i], bins=30, color=color,
                         alpha=0.7, edgecolor="black")
            axes[i].set_title(loss_type)
            axes[i].set_xlabel("Loss")
            axes[i].set_ylabel("Frequency")
            axes[i].grid(True)
        axes[3].hist(np.sum(targets, axis=1), bins=30, color=colors[3],
                     alpha=0.7, edgecolor="black")
        axes[3].set_title("Total Loss")
        plt.tight_layout()
        plt.show()


class IndivPredDataset(PredictorDataset):

    def _get_list_file_path(self, ds):
        return os.path.join(LOG_DIR, ds, "train_indiv_DELSSOME-pfic",
                            self.phase, "sub_list.txt")

    def _get_param_folder(self, ds, group_or_sub_id):
        return get_run_dir(ds, "indiv_Euler-pfic", "train", "1", "1",
                           group_or_sub_id)


class IndivPredDataModule(DataModule):

    @property
    def ds_class(self):
        return IndivPredDataset


# --------------------------------------------------------------------------- #
# Multimodal (transformer) predictor
#   inputs: (ds_idx, params[N,3], sc_euler[N,N], emp_fc[N,N], fcd_dl[bins], G)
# --------------------------------------------------------------------------- #
class MulModPredDs(PredictorDataset):

    @property
    def cache_path_prefix(self):
        return f"multimodal-{super().cache_path_prefix}"

    def param_dict_to_data(self, ds_idx, param_dict, used_n_children, emp_stats):
        data, targets = [], []
        sc_euler = emp_stats["sc_euler"]
        emp_fc = emp_stats["emp_fc"]
        fcd_dl = emp_stats["fcd_dl"]
        N = sc_euler.shape[0]
        for i, child_idx in enumerate(
                param_dict["valid_param_indices"][:used_n_children]):
            param_vector = param_dict["parameter"][:, child_idx].to(self.default_dtype)
            # Reshape the (3N+1) parameter vector into a (N, 3) per-ROI tensor
            # plus the scalar global coupling G.
            param_tensor = torch.stack(
                [param_vector[:N], param_vector[N:2 * N], param_vector[2 * N + 1:]],
                dim=1,
            )
            scale_factor = param_vector[2 * N]
            data.append((ds_idx, param_tensor, sc_euler, emp_fc, fcd_dl, scale_factor))
            corr_loss = param_dict["corr_loss"][i]
            l1_loss = param_dict["l1_loss"][i]
            ks_loss = param_dict["ks_loss"][i]
            targets.append(torch.tensor([corr_loss, l1_loss, ks_loss],
                                        dtype=self.default_dtype))
        return data, targets


class IndivMulModPredDs2014(MulModPredDs):
    """Individual-level, MFM-2014, multimodal-predictor dataset."""

    def _get_list_file_path(self, ds):
        return os.path.join(LOG_DIR, ds, "train_indiv_DELSSOME-pfic",
                            self.phase, "sub_list.txt")

    def _get_param_folder(self, ds, group_or_sub_id):
        # NOTE: assumes the Euler CMA-ES training data was generated with
        # ``trial_idx=1`` and ``seed_idx=1`` (the convention used in the paper).
        return get_run_dir(ds, "indiv_Euler-pfic", "train", "1", "1",
                           group_or_sub_id)

    @property
    def cache_path_prefix(self):
        return f"pfic-indiv-{super().cache_path_prefix}"


class IndivMulModPredDm2014(DataModule):

    @property
    def ds_class(self):
        return IndivMulModPredDs2014
