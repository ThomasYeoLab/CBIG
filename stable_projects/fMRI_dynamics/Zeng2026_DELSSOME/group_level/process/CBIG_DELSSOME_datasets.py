# Written by Tianchu Zeng and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

import torch
import torch.utils.data as Data


class DataSetPredictor(Data.Dataset):

    def __init__(self, parameter_sets, loss_sets, sc_sets, fc_sets, fcd_sets):
        """Get preprocessed tensors to generate dataset for MFM deep learning model

        Args:
            parameter_sets (tensor): [n_set, n_param, 3 * n_roi + 1]
            loss_sets (tensor): [n_set, n_param, 3]
            sc_sets (tensor): [n_set, n_roi, n_roi]
            fc_sets (tensor): [n_set, n_roi, n_roi]
            fcd_sets (tensor): [n_set, bins]
        """
        super().__init__()
        self.n_set = sc_sets.shape[0]
        n_roi = sc_sets.shape[1]
        self.n_param = parameter_sets.shape[1]  # n_param = n_epoch * (param sets per epoch)

        ut_mask = torch.ones(n_roi, n_roi, dtype=torch.bool)
        ut_mask = torch.triu(ut_mask, 1)

        sc_sets_1d = sc_sets[:, ut_mask]
        # Normalize SC to a comparable scale: 0.02 is an empirical scaling factor applied before multiplying by G.
        # The upper-triangle max is used for normalization but the full matrix is kept as the feature.
        self.sc_sets = 0.02 * sc_sets / torch.max(sc_sets_1d, dim=1, keepdim=True)[0].unsqueeze(2)

        self.fc_sets = fc_sets  # [n_param, n_roi, n_roi]
        # FCD is stored as a CDF; torch.diff with a zero prepend converts it to a PDF for use as a density input.
        tmp_prepend = torch.tensor([0]).unsqueeze(1).expand(fcd_sets.shape[0], -1)
        self.fcd_pdf_sets = torch.diff(fcd_sets, dim=1, prepend=tmp_prepend)

        self.loss_sets = loss_sets  # [n_set, n_param, 3]

        wEE = parameter_sets[:, :, 0:n_roi]  # [n_set, n_param, n_roi]
        wEI = parameter_sets[:, :, n_roi:2 * n_roi]
        G = parameter_sets[:, :, 2 * n_roi:2 * n_roi + 1]  # [n_set, n_param, 1]
        sigma = parameter_sets[:, :, 2 * n_roi + 1:]

        # concatenate wEE, wEI, sigma as [n_set, n_param, n_roi, 3]
        self.parameter_sets = torch.stack((wEE, wEI, sigma), dim=-1)  # [n_set, n_param, n_roi, 3]

        # expand the G to match the shape of sc_sets
        G = G.unsqueeze(2)
        # expand the shape of sc
        self.sc_sets = self.sc_sets.unsqueeze(1).expand(-1, G.shape[1], -1, -1)  # [n_set, n_param, n_roi, n_roi]
        self.sc_sets = G * self.sc_sets  # [n_set, n_param, n_roi, n_roi]

    def __getitem__(self, item):
        set_idx = item // self.n_param
        param_idx = item % self.n_param
        return (self.parameter_sets[set_idx, param_idx], self.loss_sets[set_idx, param_idx],
                self.sc_sets[set_idx, param_idx], self.fc_sets[set_idx], self.fcd_pdf_sets[set_idx])

    def __len__(self):
        return self.n_set * self.n_param
