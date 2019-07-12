#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Tong He and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import numpy as np
import scipy.io as sio

from cbig.He2019.config import config
from cbig.He2019.CBIG_prepare_data import data_ukbb_fnn, data_ukbb_fnn_sex
from cbig.He2019.CBIG_prepare_data import data_hcp_fnn, data_ukbb_brainnetcnn
from cbig.He2019.CBIG_prepare_data import data_ukbb_gcnn, data_ukbb_gcnn_sex
from cbig.He2019.CBIG_prepare_data import data_hcp_brainnetcnn, data_hcp_gcnn
from cbig.He2019.CBIG_prepare_data import get_gcnn_graph
from cbig.He2019.CBIG_prepare_data import data_ukbb_brainnetcnn_sex


def main():
    """main function to load and process dataset

    Returns:
        None
    """
    # load and process for HCP dataset
    os.makedirs(config.HCP_INTER_DIR, exist_ok=True)
    data_hcp_fnn(config.HCP_INTER_DIR, config.HCP_NUM_FOLD)
    data_hcp_brainnetcnn(config.HCP_INTER_DIR, config.HCP_NUM_FOLD)
    data_hcp_gcnn(config.HCP_INTER_DIR, config.HCP_NUM_FOLD)

    # load and process for UK biobank dataset
    os.makedirs(config.UKBB_INTER_DIR, exist_ok=True)
    data_ukbb_fnn(config.UKBB_INTER_DIR)
    data_ukbb_fnn_sex(config.UKBB_INTER_DIR)
    data_ukbb_brainnetcnn(config.UKBB_INTER_DIR)
    data_ukbb_brainnetcnn_sex(config.UKBB_INTER_DIR)
    data_ukbb_gcnn(config.UKBB_INTER_DIR)
    data_ukbb_gcnn_sex(config.UKBB_INTER_DIR)

    # generate adjacency matrix for GCNN
    graph_dir = config.GRAPH_FOLDER
    mat_content = sio.loadmat(
        os.path.join(config.UKBB_ORIG_DIR, config.UKBB_CORR_MAT))
    x = np.transpose(mat_content['corr_mat'], (2, 0, 1))
    get_gcnn_graph(graph_dir, x)
    get_gcnn_graph(graph_dir, x, k=1)

    mat_content = sio.loadmat(
        os.path.join(config.HCP_ORIG_DIR, config.HCP_CORR_MAT))
    x = np.transpose(mat_content['corr_mat'], (2, 0, 1))
    get_gcnn_graph(graph_dir, x)

    return


if __name__ == '__main__':
    main()
