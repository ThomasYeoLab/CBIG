# Written by Fang Tian, Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

"""
DELSSOME cost-predictor neural-network architectures.

The model that DELSSOME CMA-ES uses in this release is
:class:`MulModLossPredictorMfm` — the multimodal (params + SC + FC + FCD-PDF)
transformer-based predictor described in Fig. 9 of Zeng et al., 2026.

The class hierarchy is:

    LossPredictor (LightningModule, abstract)
        └─ LossPredictorMfm (defines MLP-only "fc/sc/fcd/param" embeddings,
                             still used by some baselines/tests)
              └─ MulModLossPredictorMfm   (transformer architecture; Fig. 9)

We keep the older MLP-only :class:`LossPredictorMfm` because it is the natural
fall-back when one wants to retrain without transformers, and it shares the
training/evaluation plumbing with the transformer model.
"""

import os
import numpy as np
import pandas as pd
import torch
import optuna
from torch import nn, optim
from lightning import LightningModule
from torchmetrics.regression import PearsonCorrCoef

from DELSSOME_indiv.CBIG_constants import ALL_TRAIN_DS


# --------------------------------------------------------------------------- #
# Utilities
# --------------------------------------------------------------------------- #
def get_log_hparams(hparams: dict) -> dict:
    """Cast list/dict-valued hyper-parameters to strings for TensorBoard."""
    log_hparams = {}
    for key, value in hparams.items():
        if isinstance(value, (list, dict)):
            log_hparams[key] = str(value)
        else:
            log_hparams[key] = value
    return log_hparams


class MLP(nn.Module):
    """Vanilla fully-connected MLP with GELU activations between layers."""

    def __init__(self, input_dim, hidden_dims, output_dim, activation=nn.GELU):
        super().__init__()
        dims = [input_dim] + hidden_dims + [output_dim]
        layers = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < len(dims) - 2:
                layers.append(activation())
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


# --------------------------------------------------------------------------- #
# Per-modality embedders (Fig. 9 of the paper)
# --------------------------------------------------------------------------- #
class ParamEmbedder(nn.Module):
    """Embed FIC parameters per region.  ``(B, R, d_param) -> (B, R, d//2)``."""

    def __init__(self, d_param, d_half, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        self.mlp = MLP(d_param, hidden_dims, d_half)
        self.norm = nn.LayerNorm(d_half)

    def forward(self, params):
        return self.norm(self.mlp(params))


class SCEmbedder(nn.Module):
    """Embed structural connectivity per region. ``(B, R, R) -> (B, R, d//2)``.

    The MLP treats every region's SC row as a feature vector (after the row is
    pre-scaled by the global coupling ``G`` inside :func:`process_batch`).
    """

    def __init__(self, n_regions, d_half, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        self.mlp = MLP(n_regions, hidden_dims, d_half)
        self.norm = nn.LayerNorm(d_half)

    def forward(self, sc):
        return self.norm(self.mlp(sc))


class FCEmbedder(nn.Module):
    """Embed empirical FC by flattening its strict upper triangle.

    ``(B, R, R) -> (B, d)``.
    """

    def __init__(self, n_regions, d, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        n_ut = n_regions * (n_regions - 1) // 2
        self.mlp = MLP(n_ut, hidden_dims, d)
        self.norm = nn.LayerNorm(d)

    def forward(self, fc):
        _, R, _ = fc.shape
        triu_indices = torch.triu_indices(R, R, offset=1)
        x = fc[:, triu_indices[0], triu_indices[1]]
        return self.norm(self.mlp(x))


class FCDPDFEmbedder(nn.Module):
    """Embed the empirical FCD-PDF. ``(B, bins) -> (B, d)``."""

    def __init__(self, bins, d, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        self.scale_norm = nn.LayerNorm(bins)
        self.mlp = MLP(bins, hidden_dims, d)
        self.norm = nn.LayerNorm(d)

    def forward(self, fcd_pdf):
        x = self.scale_norm(fcd_pdf)
        return self.norm(self.mlp(x))


# --------------------------------------------------------------------------- #
# Base LightningModule for FC+FCD cost predictors
# --------------------------------------------------------------------------- #
class LossPredictor(LightningModule):
    """Abstract base — subclasses implement ``_init_hparams`` and ``_init_models``."""

    def __init__(self,
                 optuna_trial: optuna.Trial | None = None,
                 model_hparams: dict | None = None) -> None:
        super().__init__()
        if model_hparams is None:
            model_hparams = {}
        print("Initializing LossPredictor with optuna_trial:", optuna_trial,
              "and model_hparams:", model_hparams)
        self._init_hparams(optuna_trial, model_hparams)
        # ``save_hyperparameters`` persists the dict to checkpoints so we can
        # rebuild the architecture verbatim when loading. ``optuna_trial`` is
        # excluded to avoid storing the entire study object.
        self.save_hyperparameters(ignore=["optuna_trial"])
        print("Model hyperparameters:", self.hparams.model_hparams)
        self._init_models()
        self.loss_fn = nn.MSELoss()
        self._init_metrics()

    def _init_hparams(self, optuna_trial, model_hparams=None):
        raise NotImplementedError

    def _init_models(self):
        raise NotImplementedError

    def _init_metrics(self):
        pass

    def on_fit_start(self):
        self.logger.log_hyperparams(get_log_hparams(self.hparams.model_hparams))


# --------------------------------------------------------------------------- #
# MLP-only cost predictor (kept for backwards compatibility / baselines)
# --------------------------------------------------------------------------- #
class LossPredictorMfm(LossPredictor):
    """MLP-only DELSSOME cost predictor (param + SC + FC + FCD-PDF -> 3 costs).

    Used as an Optuna-friendly baseline; the public release uses the
    transformer-based :class:`MulModLossPredictorMfm` (overrides ``_init_models``
    and ``_init_hparams``).
    """

    def _init_hparams(self, optuna_trial, model_hparams=None):
        # Optimisation
        model_hparams["lr"] = model_hparams.get(
            "lr",
            1e-3 if optuna_trial is None
            else optuna_trial.suggest_float("lr", 1e-5, 1e-1, log=True),
        )
        model_hparams["exp_lr_gamma"] = model_hparams.get(
            "exp_lr_gamma",
            0.98 if optuna_trial is None
            else optuna_trial.suggest_float("exp_lr_gamma", 0.9, 0.99),
        )

        # IO dims (not tuned, but configurable so the same code works on
        # different parcellations / FCD bin counts).
        model_hparams.setdefault("param_input_dim", 205)
        model_hparams.setdefault("sc_input_dim", 2278)
        model_hparams.setdefault("fc_input_dim", 2278)
        model_hparams.setdefault("fcd_input_dim", 10000)
        model_hparams.setdefault("fc_score_out_dim", 2)
        model_hparams.setdefault("fcd_score_out_dim", 1)

        # Architecture
        model_hparams["last_dim"] = model_hparams.get(
            "last_dim",
            111 if optuna_trial is None
            else optuna_trial.suggest_int("last_dim", 64, 256),
        )

        # Per-modality MLP depths and widths.
        for prefix, default_dim, lower, upper in [
            ("param", 88, 32, 256),
            ("sc", 984, 128, 2048),
            ("fc", 899, 128, 2048),
            ("fcd", 5835, 128, 8192),
        ]:
            n_key = f"{prefix}_n_layers"
            dims_key = f"{prefix}_hidden_dims"
            model_hparams[n_key] = model_hparams.get(
                n_key,
                1 if optuna_trial is None
                else optuna_trial.suggest_int(n_key, 1, 4),
            )
            hidden_dims = []
            for i in range(model_hparams[n_key]):
                if optuna_trial is None:
                    hidden_dims.append(default_dim)
                else:
                    hidden_dims.append(
                        optuna_trial.suggest_int(
                            f"{prefix}_hidden_dim_layer_{i}", lower, upper
                        )
                    )
            model_hparams[dims_key] = model_hparams.get(dims_key, hidden_dims)

    def _build_branch(self, in_features, n_layers, hidden_dims, last_dim,
                      use_layer_norm):
        layers = []
        for i in range(n_layers):
            out_features = last_dim if i == n_layers - 1 else hidden_dims[i]
            layers.append(nn.Linear(in_features, out_features))
            layers.append(nn.ReLU())
            layers.append(
                nn.LayerNorm(out_features) if use_layer_norm
                else nn.BatchNorm1d(out_features)
            )
            in_features = out_features
        return nn.Sequential(*layers)

    def _init_models(self):
        h = self.hparams.model_hparams
        last_dim = h["last_dim"]
        use_ln = h.get("use_layer_norm", True)
        self.param_emb = self._build_branch(
            h["param_input_dim"], h["param_n_layers"], h["param_hidden_dims"],
            last_dim, use_ln,
        )
        self.sc_emb = self._build_branch(
            h["sc_input_dim"], h["sc_n_layers"], h["sc_hidden_dims"],
            last_dim, use_ln,
        )
        self.fc_emb = self._build_branch(
            h["fc_input_dim"], h["fc_n_layers"], h["fc_hidden_dims"],
            last_dim, use_ln,
        )
        self.fcd_emb = self._build_branch(
            h["fcd_input_dim"], h["fcd_n_layers"], h["fcd_hidden_dims"],
            last_dim, use_ln,
        )
        self.output_fc_score = nn.Sequential(
            nn.Linear(last_dim, h["fc_score_out_dim"]), nn.Sigmoid())
        self.output_fcd_score = nn.Sequential(
            nn.Linear(last_dim, h["fcd_score_out_dim"]), nn.Sigmoid())

    def _init_metrics(self):
        self.corr_loss_corr = PearsonCorrCoef()
        self.l1_loss_corr = PearsonCorrCoef()
        self.ks_loss_corr = PearsonCorrCoef()

    def forward(self, param_vectors, sc_dl, fc_dl, fcd_dl):
        param_emb = self.param_emb(param_vectors)
        sc_emb = self.sc_emb(sc_dl)
        input_emb = param_emb + sc_emb
        fc_emb = self.fc_emb(fc_dl)
        fcd_emb = self.fcd_emb(fcd_dl)
        corr_l1 = self.output_fc_score(input_emb + fc_emb)
        ks = self.output_fcd_score(input_emb + fcd_emb)
        return torch.cat((corr_l1, ks), dim=1)

    def process_batch(self, batch):
        (ds_idx, param_vector, sc_dl, fc_dl, fcd_dl), y = batch
        y_hat = self(param_vector, sc_dl, fc_dl, fcd_dl)
        return {"ds_idx": ds_idx, "y_hat": y_hat, "y": y}

    def training_step(self, batch, batch_idx):
        res = self.process_batch(batch)
        loss = self.loss_fn(res["y_hat"], res["y"])
        self.log("train_loss/total_loss", loss)
        return loss

    def validation_step(self, batch, batch_idx):
        res = self.process_batch(batch)
        loss = self.loss_fn(res["y_hat"], res["y"])
        self.log("val_loss/total_loss", loss)

    def on_test_start(self):
        self.all_y_hat = []
        self.all_y = []
        self.all_ds_indices = []

    def test_step(self, batch, batch_idx):
        res = self.process_batch(batch)
        ds_idx, y_hat, y = res["ds_idx"], res["y_hat"], res["y"]
        loss = self.loss_fn(y_hat, y)
        self.log("test_loss/total_loss", loss)
        self.all_y_hat.append(y_hat.cpu().numpy())
        self.all_y.append(y.cpu().numpy())
        self.all_ds_indices.append(ds_idx.cpu().numpy())
        self.corr_loss_corr.update(y_hat[:, 0], y[:, 0])
        self.l1_loss_corr.update(y_hat[:, 1], y[:, 1])
        self.ks_loss_corr.update(y_hat[:, 2], y[:, 2])
        self.log("test_metric/corr_loss_corr", self.corr_loss_corr)
        self.log("test_metric/l1_loss_corr", self.l1_loss_corr)
        self.log("test_metric/ks_loss_corr", self.ks_loss_corr)

    def on_test_end(self):
        y_hat = np.concatenate(self.all_y_hat, axis=0)
        y = np.concatenate(self.all_y, axis=0)
        ds_indices = np.concatenate(self.all_ds_indices, axis=0)
        # Dump (ground-truth, prediction) tables for diagnostic plots.
        df = pd.DataFrame({
            "ds_name": [ALL_TRAIN_DS[idx] if idx < len(ALL_TRAIN_DS) else str(idx)
                        for idx in ds_indices],
            "ds_indices": ds_indices,
            "total_loss": np.sum(y, axis=1),
            "corr_loss": y[:, 0],
            "l1_loss": y[:, 1],
            "ks_loss": y[:, 2],
            "pred_total_loss": np.sum(y_hat, axis=1),
            "pred_corr_loss": y_hat[:, 0],
            "pred_l1_loss": y_hat[:, 1],
            "pred_ks_loss": y_hat[:, 2],
        })
        csv_path = os.path.join(self.logger.save_dir, "test_pred.csv")
        df.to_csv(csv_path, index=False)
        print(f"Test predictions saved to: {csv_path}")

    def configure_optimizers(self):
        h = self.hparams.model_hparams
        optimizer = optim.Adam(self.parameters(), lr=h["lr"])
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer, gamma=h["exp_lr_gamma"])
        return [optimizer], [scheduler]


# --------------------------------------------------------------------------- #
# Transformer cost predictor (default DELSSOME model — Fig. 9 of the paper)
# --------------------------------------------------------------------------- #
class MulModLossPredictorMfm(LossPredictorMfm):
    """
    Multimodal DELSSOME predictor (Fig. 9 of Zeng et al., 2026).

    For each ROI we concatenate a parameter embedding and a (G-scaled) SC
    embedding into a length-``d`` token; prepending a learned CLS token gives
    a sequence of length 69 (for DK68). A transformer encoder produces a
    pooled FIC-model embedding, which is additively fused with separate FC
    and FCD-PDF embeddings to predict the three FC+FCD cost components.
    """

    def _init_hparams(self, optuna_trial, model_hparams=None):
        # Optimisation
        model_hparams["lr"] = model_hparams.get(
            "lr",
            1e-3 if optuna_trial is None
            else optuna_trial.suggest_float("lr", 1e-5, 1e-1, log=True),
        )
        model_hparams["exp_lr_gamma"] = model_hparams.get(
            "exp_lr_gamma",
            0.98 if optuna_trial is None
            else optuna_trial.suggest_float("exp_lr_gamma", 0.9, 0.99),
        )

        # Architecture (defaults match the released checkpoint)
        model_hparams.setdefault("n_regions", 68)
        model_hparams.setdefault("d_param", 3)        # wEE, wEI, sigma per ROI
        model_hparams.setdefault("bins", 10000)
        model_hparams.setdefault("fc_score_out_dim", 2)
        model_hparams.setdefault("fcd_score_out_dim", 1)

        model_hparams["d"] = model_hparams.get(
            "d",
            64 if optuna_trial is None
            else optuna_trial.suggest_categorical("d", [32, 64, 128, 256]),
        )
        model_hparams["transformer_layers"] = model_hparams.get(
            "transformer_layers",
            4 if optuna_trial is None
            else optuna_trial.suggest_int("transformer_layers", 2, 8),
        )
        model_hparams["transformer_heads"] = model_hparams.get(
            "transformer_heads",
            8 if optuna_trial is None
            else optuna_trial.suggest_categorical("transformer_heads", [4, 8, 16]),
        )
        model_hparams["mlp_n_layers"] = model_hparams.get(
            "mlp_n_layers",
            1 if optuna_trial is None
            else optuna_trial.suggest_int("mlp_n_layers", 0, 3),
        )

        mlp_hidden = []
        for i in range(model_hparams["mlp_n_layers"]):
            if optuna_trial is None:
                mlp_hidden.append(64)
            else:
                mlp_hidden.append(
                    optuna_trial.suggest_int(f"mlp_hidden_dim_layer_{i}", 32, 256)
                )
        model_hparams["mlp_hidden"] = model_hparams.get("mlp_hidden", mlp_hidden)

        # Make sure d and heads are compatible.
        if model_hparams["d"] % 2 != 0:
            model_hparams["d"] += 1
        while model_hparams["d"] % model_hparams["transformer_heads"] != 0:
            if optuna_trial is None:
                model_hparams["d"] += 1
            else:
                valid_heads = [h for h in [4, 8, 16]
                               if model_hparams["d"] % h == 0]
                if valid_heads:
                    model_hparams["transformer_heads"] = valid_heads[0]
                else:
                    model_hparams["d"] += 1

    def _init_models(self):
        h = self.hparams.model_hparams
        n_regions = h["n_regions"]
        d_param = h["d_param"]
        d = h["d"]
        bins = h["bins"]
        mlp_hidden = h["mlp_hidden"]
        assert d % 2 == 0, "Embedding dim d must be even"
        d_half = d // 2

        # Per-modality embedders.
        self.param_embed = ParamEmbedder(d_param, d_half, mlp_hidden)
        self.sc_embed = SCEmbedder(n_regions, d_half, mlp_hidden)
        self.fc_embed = FCEmbedder(n_regions, d, mlp_hidden)
        self.fcd_pdf_embed = FCDPDFEmbedder(bins, d, mlp_hidden)

        # Learned [CLS] token + transformer encoder.
        self.cls_token = nn.Parameter(torch.randn(1, 1, d))
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d,
            nhead=h["transformer_heads"],
            dim_feedforward=d * 4,
            activation="gelu",
            batch_first=True,
        )
        self.transformer = nn.TransformerEncoder(
            encoder_layer, num_layers=h["transformer_layers"])
        self.transformer_norm = nn.LayerNorm(d)

        # Fusion heads -> 2 FC scores (corr_loss, l1_loss) + 1 FCD score (ks_loss).
        self.fuse_norm_fc = nn.LayerNorm(d)
        self.fuse_mlp_fc = MLP(d, mlp_hidden, 2)
        self.fuse_norm_fcd_pdf = nn.LayerNorm(d)
        self.fuse_mlp_fcd_pdf = MLP(d, mlp_hidden, 1)
        self.output_layer = nn.Sigmoid()

    def process_batch(self, batch):
        (ds_idx, params, sc_euler, emp_fc, emp_fcd_pdf, G), y = batch
        # Scale per-region SC profiles by the (sampled) global coupling G,
        # so the SC tokens reflect the effective network strength.
        sc = sc_euler * (G.unsqueeze(-1).unsqueeze(-1))  # [B, R, R]
        y_hat = self(params, sc, emp_fc, emp_fcd_pdf)
        return {"ds_idx": ds_idx, "y_hat": y_hat, "y": y}

    def forward(self, params, sc, emp_fc, emp_fcd_pdf):
        """params: [B, R, d_param];  sc: [B, R, R];
           emp_fc: [B, R, R];  emp_fcd_pdf: [B, bins].
           Returns predicted FC+FCD costs: [B, 3]."""
        B, R, _ = params.shape
        p_emb = self.param_embed(params)             # [B, R, d/2]
        s_emb = self.sc_embed(sc)                    # [B, R, d/2]
        nodes = torch.cat([p_emb, s_emb], dim=-1)    # [B, R, d]
        cls_tokens = self.cls_token.expand(B, -1, -1)
        seq = torch.cat([cls_tokens, nodes], dim=1)  # [B, R+1, d]
        t_out = self.transformer(seq)
        t_out = self.transformer_norm(t_out)
        cls_out = t_out[:, 0, :]                     # [B, d]

        fc_out = self.fc_embed(emp_fc)
        joint_fc = self.fuse_norm_fc(cls_out + fc_out)
        out_fc = self.fuse_mlp_fc(joint_fc)          # [B, 2]

        fcd_pdf_out = self.fcd_pdf_embed(emp_fcd_pdf)
        joint_fcd_pdf = self.fuse_norm_fcd_pdf(cls_out + fcd_pdf_out)
        out_fcd_pdf = self.fuse_mlp_fcd_pdf(joint_fcd_pdf)  # [B, 1]

        # Sigmoid-bound costs to [0, 1].
        out_fc = self.output_layer(out_fc)
        out_fcd_pdf = self.output_layer(out_fcd_pdf)
        return torch.cat([out_fc, out_fcd_pdf], dim=-1)
