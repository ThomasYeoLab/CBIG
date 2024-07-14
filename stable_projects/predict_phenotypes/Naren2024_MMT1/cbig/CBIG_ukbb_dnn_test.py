#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import argparse
import random
import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from cbig.CBIG_model_pytorch import vol_dataset
from cbig.CBIG_mics import read_datapath, load_data
from cbig.config import config


def test(args):
    '''main function for testing base DNN network (used for
       generating input for meta-matching stacking)

    Args:
        args: args from command line

    Returns:
        None

    '''

    # set all the seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.enabled = True

    # set gpu
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # load data
    icv_test, tes_sub_list, tes_phe_name, y_test = load_data(args)
    x_test = read_datapath(args.data_dir, list(map(str, tes_sub_list[0])),
                           args.dataset)
    n_phe = args.outputdim

    dset_test = vol_dataset(x_test, y_test, icv=icv_test)
    batch_size = args.batch_size
    testLoader = DataLoader(dset_test,
                            batch_size=batch_size,
                            shuffle=False,
                            num_workers=batch_size)

    # load trained model
    opt_index = args.index
    weight_path = os.path.join(
        args.model_dir, 'dnn_model_save_base',
        'CBIG_ukbb_dnn_run_0_epoch_' + str(opt_index) + '.pkl_torch')

    net = torch.load(weight_path)  # map_location=torch.device('cpu')
    net.to(device)
    net.train(False)

    record_pred = np.zeros((0, n_phe))  # prediction value
    tes_res_record = np.zeros((1, 1, x_test.shape[0], n_phe))
    for (x, y, icv) in testLoader:
        x, y, icv = x.to(device), y.to(device), icv.to(device)
        outputs = net(x, icv)
        record_pred = np.concatenate((record_pred, outputs.data.cpu().numpy()),
                                     axis=0)

        del outputs, x, y
        torch.cuda.empty_cache()
    tes_res_record[0, 0, :, :] = np.squeeze(record_pred)

    log_args = {'tes_res_record': tes_res_record}

    # save record value for future use
    file_str = 'dnn_test_base.npz'
    name_str = os.path.join(args.out_dir, args.out_subdir, file_str)
    os.makedirs(os.path.join(args.out_dir, args.out_subdir), exist_ok=True)
    np.savez(name_str, **log_args)

    return


def get_args():
    '''function to get args from the command line and return the args

    Returns:
        argparse.ArgumentParser: args that could be used by other function
    '''
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--phe_dir', type=str, default=None)
    parser.add_argument('--icv_dir', type=str, default=None)
    parser.add_argument('--sub_dir', type=str, default=None)
    parser.add_argument('--tra_sub_dir', type=str, default=None)
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--model_dir', type=str, default=None)
    parser.add_argument('--dataset', type=str, default=config.DATASET)
    parser.add_argument('--out_dir', '-o', type=str, default=None)
    parser.add_argument('--out_subdir', '-osub', type=str, default=None)
    parser.add_argument('--seed', type=int, default=config.SEED)
    parser.add_argument('--outputdim', type=int, default=None)
    parser.add_argument('--batch_size', type=int, default=None)
    parser.add_argument("--ukbb_icv_dir", type=str, default=None)
    parser.add_argument('--index', type=int, default=None)
    parser.add_argument("--start_idx", type=int, default=None)
    parser.add_argument("--end_idx", type=int, default=None)
    across_dataset_parser = parser.add_mutually_exclusive_group(required=False)
    across_dataset_parser.add_argument('--across-dataset',
                                       dest='across_dataset',
                                       action='store_true')
    across_dataset_parser.add_argument('--not-across-dataset',
                                       dest='across_dataset',
                                       action='store_false')
    parser.set_defaults(across_dataset=False)

    return parser.parse_args()


if __name__ == '__main__':
    test(get_args())
