"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import torch
import torch.nn as nn


class MLP(nn.Module):
    """
    Simple fully-connected MLP with configurable layers.
    """

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


class ParamEmbedder(nn.Module):
    """
    Embed biophysical parameters per region via an MLP bottleneck.
    Input: (B, R, d_param) -> Output: (B, R, d//2)
    """

    def __init__(self, d_param, d_half, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        self.mlp = MLP(d_param, hidden_dims, d_half)
        self.norm = nn.LayerNorm(d_half)

    def forward(self, params):
        # params: [B, R, d_param]
        B, R, _ = params.shape
        x = self.mlp(params)  # [B, R, d_half]
        x = self.norm(x)
        return x


class SCEmbedder(nn.Module):
    """
    Embed structural connectivity per region via MLP or GNN bottleneck.
    Input: (B, R, R) -> Output: (B, R, d//2)
    """

    def __init__(self, n_regions, d_half, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        # treat each row of SC adjacency as features
        self.mlp = MLP(n_regions, hidden_dims, d_half)
        self.norm = nn.LayerNorm(d_half)

    def forward(self, sc):
        # sc: [B, R, R]
        B, R, _ = sc.shape
        x = self.mlp(sc)  # [B, R, d_half]
        x = self.norm(x)
        return x


class FCEmbedder(nn.Module):
    """
    Embed empirical FC matrix by flattening upper triangle and an MLP.
    Input: (B, R, R) -> Output: (B, d)
    """

    def __init__(self, n_regions, d, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        # number of upper-triangular entries
        n_ut = n_regions * (n_regions - 1) // 2
        self.mlp = MLP(n_ut, hidden_dims, d)
        self.norm = nn.LayerNorm(d)

    def forward(self, fc):
        # fc: [B, R, R]
        B, R, _ = fc.shape
        # extract upper triangle indices
        triu_indices = torch.triu_indices(R, R, offset=1)
        x = fc[:, triu_indices[0], triu_indices[1]]  # [B, n_ut]
        x = self.mlp(x)  # [B, d]
        x = self.norm(x)
        return x


class FCDPDFEmbedder(nn.Module):
    """
    Embed empirical FCD matrix PDF
    Input: (B, bins) -> Output: (B, d)
    """

    def __init__(self, bins, d, hidden_dims=None):
        super().__init__()
        hidden_dims = hidden_dims or []
        self.scale_norm = nn.LayerNorm(bins)
        self.mlp = MLP(bins, hidden_dims, d)
        self.norm = nn.LayerNorm(d)

    def forward(self, fcd_pdf):
        # fcd_pdf: [B, bins]
        x = self.scale_norm(fcd_pdf)  # Normalize bins
        x = self.mlp(x)  # [B, d]
        x = self.norm(x)
        return x


class MultiModalCostPredictor(nn.Module):
    """
    Multimodal DELSSOME cost predictor:
      - Embeds parameters and SC per region, fuses via Transformer
      - Embeds empirical FC directly
      - Fuses both streams to predict a scalar correlation
    """

    def __init__(
        self,
        n_regions: int,
        d_param: int,
        d: int,
        bins: int = 10000,
        transformer_layers: int = 4,
        transformer_heads: int = 8,
        mlp_hidden: list = None,
    ):
        super().__init__()
        assert d % 2 == 0, "Embedding dim d must be even"
        d_half = d // 2
        mlp_hidden = mlp_hidden or []

        # modality embedders
        self.param_embed = ParamEmbedder(d_param, d_half, mlp_hidden)
        self.sc_embed = SCEmbedder(n_regions, d_half, mlp_hidden)
        self.fc_embed = FCEmbedder(n_regions, d, mlp_hidden)
        self.fcd_pdf_embed = FCDPDFEmbedder(bins, d, mlp_hidden)

        # CLS token
        self.cls_token = nn.Parameter(torch.randn(1, 1, d))

        # Transformer Encoder
        encoder_layer = nn.TransformerEncoderLayer(d_model=d,
                                                   nhead=transformer_heads,
                                                   dim_feedforward=d * 4,
                                                   activation='gelu',
                                                   batch_first=True)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=transformer_layers)
        self.transformer_norm = nn.LayerNorm(d)

        # fusion MLP to scalar
        self.fuse_norm_fc = nn.LayerNorm(d)
        self.fuse_mlp_fc = MLP(d, mlp_hidden, 2)

        self.fuse_norm_fcd_pdf = nn.LayerNorm(d)
        self.fuse_mlp_fcd_pdf = MLP(d, mlp_hidden, 1)

        self.output_layer = nn.Sigmoid()  # Ensure output is between 0 and 1

    def forward(self, params, sc, emp_fc, emp_fcd_pdf):
        """
        params: [B, R, d_param]
        sc:     [B, R, R]
        emp_fc: [B, R, R]
        returns: [B, 1]  predicted Pearson r
        """
        B, R, _ = params.shape
        # embed streams
        p_emb = self.param_embed(params)  # [B, R, d/2]
        s_emb = self.sc_embed(sc)  # [B, R, d/2]
        # fuse node embeddings
        nodes = torch.cat([p_emb, s_emb], dim=-1)  # [B, R, d]
        # prepend CLS
        cls_tokens = self.cls_token.expand(B, -1, -1)  # [B, 1, d]
        seq = torch.cat([cls_tokens, nodes], dim=1)  # [B, R+1, d]
        # transformer
        t_out = self.transformer(seq)  # [B, R+1, d]
        t_out = self.transformer_norm(t_out)
        cls_out = t_out[:, 0, :]  # [B, d]

        # embed empirical FC
        fc_out = self.fc_embed(emp_fc)  # [B, d]

        # fuse and predict
        joint_fc = self.fuse_norm_fc(cls_out + fc_out)  # [B, d]
        out_fc = self.fuse_mlp_fc(joint_fc)  # [B, 2]

        # embed FCD PDF
        fcd_pdf_out = self.fcd_pdf_embed(emp_fcd_pdf)  # [B, d]
        joint_fcd_pdf = self.fuse_norm_fcd_pdf(cls_out + fcd_pdf_out)
        out_fcd_pdf = self.fuse_mlp_fcd_pdf(joint_fcd_pdf)  # [B, 1]

        # For FC outputs, scale to [0, 2]
        out_fc = self.output_layer(out_fc) * 2  # Scale to [0, 2]
        # For FCD PDF outputs, scale to [0, 1]
        out_fcd_pdf = self.output_layer(out_fcd_pdf)
        # concatenate outputs
        out = torch.cat([out_fc, out_fcd_pdf], dim=-1)

        return out
