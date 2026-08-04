# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
Optuna + Lightning plumbing for training the DELSSOME cost predictor.

This is a slimmed-down port of ``utils/dl_train_utils.py`` (only the model
types used in the paper's individual-level pipeline are kept). The flow is:

1. :func:`tune` builds an Optuna study, picks hyper-parameters per trial, and
   for each trial wires up a ``LightningDataModule`` and a ``LightningModule``.
2. Each trial runs ``Trainer.fit`` (up to ``max_epochs``) and returns the best
   validation loss; Optuna then drives the search.

To save GPU/CPU on a laptop set ``n_trials=1`` and ``max_epochs`` to a small
value; the paper's settings (50 trials, 15 epochs) need a GPU and about 65 hours.
"""

import os
from pathlib import Path
import pickle
import warnings

import torch
from torch import nn  # noqa: F401 — useful in user customisations
import optuna
from optuna.samplers import TPESampler
from optuna.trial import TrialState
from optuna_integration import PyTorchLightningPruningCallback
from optuna.visualization import (
    plot_slice, plot_optimization_history, plot_param_importances,
)
from lightning import Trainer
from lightning.pytorch.callbacks import ModelCheckpoint, LearningRateMonitor
from lightning.pytorch.loggers import TensorBoardLogger

from DELSSOME_indiv.datasets.CBIG_dl_dataset import (
    IndivPredDataModule, IndivMulModPredDm2014,
)
from DELSSOME_indiv.models.CBIG_DELSSOME_predictor import (
    LossPredictorMfm, MulModLossPredictorMfm,
)


# --------------------------------------------------------------------------- #
# CLI helpers
# --------------------------------------------------------------------------- #
def get_n_sample_per_phase(n_sample_per_phase_str):
    """Parse ``"train40_val20_test20-train40_val20"`` into a list of dicts."""
    n_sample_per_phase = []
    for dict_str in n_sample_per_phase_str.split("-"):
        curr_dict = {}
        for phase_n_sample in dict_str.split("_"):
            if phase_n_sample.startswith("train"):
                curr_dict["train"] = int(phase_n_sample[5:])
            elif phase_n_sample.startswith("val"):
                curr_dict["val"] = int(phase_n_sample[3:])
            elif phase_n_sample.startswith("test"):
                curr_dict["test"] = int(phase_n_sample[4:])
        n_sample_per_phase.append(curr_dict)
    return n_sample_per_phase


def walltime_to_seconds(walltime):
    """Convert ``HH:MM:SS`` to seconds (used as the Optuna study timeout)."""
    hours, minutes, seconds = map(int, walltime.split(":"))
    return hours * 3600 + minutes * 60 + seconds


def config_gpu():
    """Print visible CUDA devices and enable TF32 matmul if a GPU is present."""
    for i in range(torch.cuda.device_count()):
        info = torch.cuda.get_device_properties(i)
        print(f"CUDA:{i} {info.name}, {info.total_memory / 1024 ** 2}MB")
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision("high")


# --------------------------------------------------------------------------- #
# Study persistence
# --------------------------------------------------------------------------- #
def get_study(save_dir, seed):
    """Resume an Optuna study from ``save_dir/study.pkl`` (or create a new one)."""
    study_file_path = os.path.join(save_dir, "study.pkl")
    if os.path.isfile(study_file_path):
        return pickle.load(open(study_file_path, "rb"))
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    sampler = TPESampler(seed=seed)
    return optuna.create_study(direction="minimize", sampler=sampler)


# --------------------------------------------------------------------------- #
# Per-trial objective
# --------------------------------------------------------------------------- #
def objective(model_cls, trial, datamodule, save_dir, max_epochs,
              is_trial_run=False, enable_progress_bar=False):
    model = model_cls(optuna_trial=trial)
    # ``torch.compile`` gives ~10-20% speed-up but requires PyTorch 2.x; comment
    # out if it causes trouble on older versions.
    model = torch.compile(model)
    train_loader = datamodule.train_dataloader()
    trainer = get_trainer(trial=trial, save_dir=save_dir, max_epochs=max_epochs,
                          is_trial_run=is_trial_run,
                          n_train_batches=len(train_loader),
                          enable_progress_bar=enable_progress_bar)
    trainer.fit(
        model,
        train_dataloaders=train_loader,
        val_dataloaders=datamodule.val_dataloader(),
    )
    return trainer.callback_metrics["val_loss/total_loss"].item()


# --------------------------------------------------------------------------- #
# Model / data-module factory
# --------------------------------------------------------------------------- #
def get_data_and_model_cls(model_type, is_trial_run=False,
                           n_sample_per_phase=None, **kwargs):
    """Return ``(DataModule_instance, LightningModule_class)`` for ``model_type``.

    Supported ``model_type`` values:

    * ``"pfic-indiv-pred-multimodal"``  — default; transformer DELSSOME
      predictor used for the paper's individual-level results.
    * ``"pfic-indiv-pred"``             — MLP-only baseline.
    """
    n_sample_per_phase = n_sample_per_phase or [{"train": 1, "val": 1, "test": 1}]
    if model_type == "pfic-indiv-pred-multimodal":
        datamodule = IndivMulModPredDm2014(
            is_trial_run=is_trial_run,
            n_sample_per_phase=n_sample_per_phase,
            **kwargs,
        )
        model_cls = MulModLossPredictorMfm
    elif model_type == "pfic-indiv-pred":
        datamodule = IndivPredDataModule(
            is_trial_run=is_trial_run,
            n_sample_per_phase=n_sample_per_phase,
            **kwargs,
        )
        model_cls = LossPredictorMfm
    else:
        raise ValueError(
            f"model_type {model_type!r} is not implemented in the public release "
            f"(use 'pfic-indiv-pred-multimodal' or 'pfic-indiv-pred')."
        )
    print(f"Using DataModule: {datamodule.__class__.__name__} "
          f"and Model: {model_cls.__name__} for model_type: {model_type}")
    return datamodule, model_cls


def model_type_to_target(model_type):
    """Map a ``model_type`` to the log-directory ``target`` used by the runner."""
    if model_type in ("pfic-indiv-pred", "pfic-indiv-pred-multimodal"):
        return "train_indiv_DELSSOME-pfic"
    raise ValueError(f"model_type {model_type!r} is not implemented.")


# --------------------------------------------------------------------------- #
# Tune + Test entry points
# --------------------------------------------------------------------------- #
def tune(model_type, max_epochs, n_trials, timeout, save_dir, seed,
         batch_size=4096, num_workers=0, is_trial_run=False,
         n_sample_per_phase=None, ds_names=None,
         enable_progress_bar=False, vis_study=True):
    """Run an Optuna study on top of Lightning trainer.fit."""
    n_sample_per_phase = n_sample_per_phase or [{"train": 1, "val": 1, "test": 1}]
    extra = {} if ds_names is None else {"ds_names": ds_names}
    datamodule, model_cls = get_data_and_model_cls(
        model_type=model_type, is_trial_run=is_trial_run,
        n_sample_per_phase=n_sample_per_phase,
        batch_size=batch_size, num_workers=num_workers,
        **extra,
    )
    datamodule.setup("fit")
    study = get_study(save_dir, seed)
    print(f"batch_size for DataModule: {datamodule.hparams.batch_size}")
    study.optimize(
        lambda trial: objective(model_cls, trial, datamodule, save_dir,
                                max_epochs, is_trial_run, enable_progress_bar),
        n_trials=n_trials,
        timeout=timeout,
        catch=[torch.cuda.OutOfMemoryError],
    )
    print_best_trial(study)
    with open(os.path.join(save_dir, "study.pkl"), "wb") as f:
        pickle.dump(study, f)
    if vis_study:
        visualize_study(study, save_dir)
    return study


def visualize_study(study, save_dir):
    """Dump Optuna's diagnostic plots to ``save_dir/*.png``.

    Plot export goes through plotly → kaleido, which since kaleido v1
    requires a headless Chrome install. On machines without Chrome the
    export raises; treat that as non-fatal — the plots are diagnostic
    and the rest of the pipeline (study.pkl, best_value.pth) is what
    downstream phases consume.
    """
    def _safe_write(fig_fn, filename):
        try:
            fig_fn().write_image(os.path.join(save_dir, filename))
        except Exception as e:
            warnings.warn(f"Skipping {filename}: {type(e).__name__}: {e}")

    _safe_write(lambda: plot_slice(study), "hparams_slice_plot.png")
    _safe_write(lambda: plot_optimization_history(study),
                "hparams_optimization_history.png")
    n_completed = sum(t.state == TrialState.COMPLETE for t in study.trials)
    if n_completed >= 2:
        _safe_write(lambda: plot_param_importances(study),
                    "hparams_param_importances.png")


def print_best_trial(study):
    print(f"Number of finished trials: {len(study.trials)}")
    trial = study.best_trial
    print("Best trial:")
    print(f"  Trial number: {trial.number}")
    print(f"  Value: {trial.value}")
    print("  Params:")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")


def get_trainer(trial, save_dir, max_epochs, is_trial_run=False,
                n_train_batches=None, **kwargs):
    """Build a :class:`lightning.Trainer` with the right callbacks/logger."""
    trial_run_args = {
        "limit_train_batches": 2,
        "limit_val_batches": 2,
        "limit_test_batches": 2,
    } if is_trial_run else {}
    # Clamp the logging cadence to the number of training batches so small
    # datasets (e.g. the released example) don't trip Lightning's
    # "log_every_n_steps is larger than the number of training batches" warning.
    if is_trial_run:
        log_every_n_steps = 1
    elif n_train_batches:
        log_every_n_steps = min(50, n_train_batches)
    else:
        log_every_n_steps = 50
    return Trainer(
        callbacks=get_trainer_callbacks(trial),
        logger=TensorBoardLogger(save_dir=save_dir),
        max_epochs=max_epochs,
        accelerator="auto",
        devices=1,
        log_every_n_steps=log_every_n_steps,
        precision="32-true",
        **trial_run_args,
        **kwargs,
    )


def get_trainer_callbacks(trial):
    """LR monitor + best-checkpoint saver + Optuna pruning."""
    callbacks = []
    callbacks.append(LearningRateMonitor(logging_interval="epoch"))
    checkpoint_callback = ModelCheckpoint(
        save_top_k=1,
        monitor="val_loss/total_loss",
        mode="min",
        save_last=True,
        filename="epoch={epoch}-step={step}-val_loss={val_loss/total_loss:.7f}",
        auto_insert_metric_name=False,
    )
    checkpoint_callback.CHECKPOINT_NAME_LAST = (
        "epoch={epoch}-step={step}-last-val_loss={val_loss/total_loss:.7f}"
    )
    callbacks.append(checkpoint_callback)
    if trial is not None:
        callbacks.append(
            PyTorchLightningPruningCallback(trial, monitor="val_loss/total_loss")
        )
    return callbacks


def test(model, model_type, save_dir,
         batch_size=4096, num_workers=0,
         is_trial_run=False,
         n_sample_per_phase=None, ds_names=None,
         enable_progress_bar=False):
    """Load a trained checkpoint (or a model instance) and run ``trainer.test``."""
    n_sample_per_phase = n_sample_per_phase or [{"train": 1, "val": 1, "test": 1}]
    extra = {} if ds_names is None else {"ds_names": ds_names}
    datamodule, model_cls = get_data_and_model_cls(
        model_type=model_type, is_trial_run=is_trial_run,
        n_sample_per_phase=n_sample_per_phase,
        batch_size=batch_size, num_workers=num_workers,
        **extra,
    )
    if isinstance(model, str):
        model = model_cls.load_from_checkpoint(model)
    trainer = get_trainer(trial=None, save_dir=save_dir, max_epochs=1,
                          is_trial_run=is_trial_run,
                          enable_progress_bar=enable_progress_bar)
    trainer.test(model, datamodule=datamodule)
    return model
