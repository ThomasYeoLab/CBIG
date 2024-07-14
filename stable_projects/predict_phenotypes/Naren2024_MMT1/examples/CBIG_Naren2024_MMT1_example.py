#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import os
import random
import argparse
import numpy as np
import pandas as pd
from sklearn.datasets import make_regression
import nibabel as nib
from nibabel.testing import data_path
import subprocess


def example_args_parser(setup):
    '''function to get args from command line and return the args

    Returns:
       argparse.ArgumentParser: args that could be used by other function
    '''

    parser = argparse.ArgumentParser(prog='ExampleArgs')
    parser.add_argument('--n_samples', type=int, default=300)
    parser.add_argument('--n_tra_samples', type=int, default=50)
    parser.add_argument('--n_dim', type=int, default=10)
    parser.add_argument('--n_shape',
                        nargs='+',
                        type=int,
                        default=[182, 218, 182])
    parser.add_argument('--n_informative', type=int, default=7)
    parser.add_argument('--n_targets', type=int, default=5)
    parser.add_argument('--setup', type=str, default=setup)
    parser.add_argument('--n_tra_targets', type=int, default=3)
    parser.add_argument('--base_dir',
                        type=str,
                        default=os.path.abspath(os.path.join(os.getcwd())))

    args, _ = parser.parse_known_args()

    return args


def save_meta_set(dir_out, eid, phe, idps, ind, x_3d, icv, df_y, df_x, name):
    '''function to save generated fake data

    Args:
        dir_out (str): save directory
        eid (ndarray): subject id
        phe (ndarray): phenotype name
        idps (ndarray) : IDP name
        ind (ndarray) : index for training or testing
        x_3d (ndarray) : volumetric data for 3D CNN Model
        icv (ndarray) : icv data
        df_y (ndarray) : output data
        df_x (ndarray) : input data for elasticnet
        name (str) : save file name

    Returns:
           None
    '''

    np.savetxt(os.path.join(dir_out, name + '_subj_list.txt'), eid, fmt='%s')
    np.savetxt(os.path.join(dir_out, name + '_final_phe_list.txt'),
               phe,
               fmt='%s')
    df_tmp = df_y.copy()
    df_tmp = df_tmp[['eid'] + phe]
    df_tmp = df_tmp.loc[df_tmp['eid'].isin(eid)]
    df_tmp.to_csv(os.path.join(dir_out, name + '_final.csv'), index=False)

    df_tmp = df_x.copy()
    df_tmp = df_tmp[['eid'] + idps]
    df_tmp = df_tmp.loc[df_tmp['eid'].isin(eid)]
    df_tmp.to_csv(os.path.join(dir_out, name + '_idps.csv'), index=False)

    df_tmp = pd.DataFrame(icv[ind], columns=['inverse_determinant'])
    df_tmp['eid'] = eid
    cols = df_tmp.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_tmp = df_tmp[cols]
    df_tmp.to_csv(os.path.join(dir_out, name + '_icv.csv'), index=False)

    example_filename = os.path.join(data_path, 'example4d.nii.gz')
    exp_img = nib.load(example_filename)

    for i in range(len(eid)):
        img = nib.Nifti1Image(x_3d[ind[i], :], affine=exp_img.affine)
        os.makedirs(os.path.join(dir_out, 'T1', str(eid[i])), exist_ok=True)
        nib.save(
            img,
            os.path.join(dir_out, 'T1', str(eid[i]),
                         "T1_MNILinear_1mm_newsize.nii.gz"))


def generate_data(args):
    '''function to generate fake data

    Args:
        args: args from command line

    Returns:
        None
    '''

    base_dir = args.base_dir
    dir_out = os.path.join(base_dir, args.setup, 'data')
    os.makedirs(dir_out, exist_ok=True)

    # parameter for data generation
    n_samples = args.n_samples
    n_tra_samples = args.n_tra_samples
    n_dim = args.n_dim  # for elasticnet
    n_shape = args.n_shape  # for 3D cnn
    n_informative = args.n_informative
    n_targets = args.n_targets
    n_tra_targets = args.n_tra_targets

    # generate data for elasticnet
    x, y = make_regression(n_samples=n_samples,
                           n_features=n_dim,
                           n_informative=n_informative,
                           n_targets=n_targets,
                           noise=0.1,
                           random_state=0)

    print('randomly generated x with shape', x.shape, 'y with shape', y.shape)

    # convert x to 3D format
    icv, x_res = x[:, 0], x[:, 1:]

    t_dim = int(n_shape[0] * n_shape[1] * n_shape[2])
    r_num = int(round(t_dim / (n_dim - 1), 0))
    x_3d = np.tile(x_res, (1, r_num))
    x_3d = x_3d[:, :t_dim].reshape(
        (n_samples, n_shape[0], n_shape[1], n_shape[2]))
    x_3d = (x_3d + 0.5) * 0.1

    print('converted to fc with shape', x_3d.shape)

    # generate eid and idp name
    eid = np.arange(n_samples) + 8000001
    idps = [str(901 + i) + '-999.0' for i in np.arange(n_dim)]

    # put x into dataframe with eid and name
    df_x = pd.DataFrame(x, columns=idps)
    df_x['eid'] = eid
    cols = df_x.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_x = df_x[cols]

    # generate eid and phenotype name
    phes = [str(101 + i) + '-888.0' for i in np.arange(n_targets)]

    # put y into dataframe with eid and name
    df_y = pd.DataFrame(y, columns=phes)
    df_y['eid'] = eid
    cols = df_y.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_y = df_y[cols]

    # split data and phenotype
    ind_tra = range(0, n_tra_samples)
    ind_tes = range(n_tra_samples, n_samples)
    eid_tra = eid[ind_tra]
    eid_tes = eid[ind_tes]
    phe_tra = phes[:n_tra_targets]
    phe_tes = phes[n_tra_targets:]

    # save out data generated
    save_meta_set(dir_out, eid_tes, phe_tes, idps, ind_tes, x_3d, icv, df_y,
                  df_x, args.setup + '_test')
    save_meta_set(dir_out, eid_tra, phe_tra, idps, ind_tra, x_3d, icv, df_y,
                  df_x, args.setup + '_train')


def test_elasticnet(args):
    '''function to test elasticnet function
    Args:
        args: args from command line

    Returns:
        None
    '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    dataset = args.setup + "_test_dataset"

    in_dir = os.path.join(args.base_dir, args.setup, 'data')
    phe_dir = os.path.join(in_dir, args.setup + '_test' + '_final.csv')
    sub_dir = os.path.join(in_dir, args.setup + '_test' + '_subj_list.txt')
    data_dir = os.path.join(in_dir, args.setup + '_test' + '_idps.csv')
    out_dir = os.path.join(res_dir, args.setup + "_test_dataset")
    out_subdir = "/elasticnet"
    inter_dir = os.path.join(out_dir, "output_intermediate")
    n_rng = "2"
    phe_start_idx = "1"
    phe_end_idx = "3"
    num_inner_folds = "5"
    ks = ["10", "20", "50", "100", "200"]
    seed = "0"

    subprocess.call([
        'python', '-m', 'utils.generate_split_kshot', '--phe_dir', phe_dir,
        '--sub_dir', sub_dir, '--inter_dir', inter_dir, '--n_rng', n_rng,
        '--start_idx', phe_start_idx, '--end_idx', phe_end_idx,
        '--num_inner_folds', num_inner_folds, '--ks', *ks, '--dataset',
        dataset, '--seed', seed
    ])

    subprocess.call([
        'python', '-m', 'ElasticNet.elasticnet', '--phe_dir', phe_dir,
        '--sub_dir', sub_dir, '--data_dir', data_dir, '--out_dir', out_dir,
        '--out_subdir', out_subdir, '--inter_dir', inter_dir, '--rng', n_rng,
        '--across-dataset', '--start_idx', phe_start_idx, '--end_idx',
        phe_end_idx, '--ks', *ks, '--dataset', dataset
    ])


def test_DNN_MM_training(args):
    '''function to test DNN training function
    Args:
        args: args from command line

    Returns:
        None
    '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    dataset = args.setup + "_train_dataset"

    in_dir = os.path.join(args.base_dir, args.setup, 'data')
    phe_dir = os.path.join(in_dir, args.setup + '_train' + '_final.csv')
    sub_dir = os.path.join(in_dir, args.setup + '_train' + '_subj_list.txt')
    icv_dir = os.path.join(in_dir, args.setup + '_train' + '_icv.csv')
    data_dir = os.path.join(in_dir, 'T1')
    out_dir = os.path.join(res_dir, args.setup + "_train_dataset")
    out_subdir = "trained_model_ukbb"
    inter_dir = os.path.join(out_dir, "output_intermediate")
    epochs = "2"
    batch_size = "50"
    lr = "1e-5"
    scheduler_decrease = "1"
    cn = ["1", "2", "2", "3", "2", "2"]
    output_dim = "3"
    dropout = "0.1"
    phe_start_idx = "1"
    phe_end_idx = "4"

    subprocess.call([
        'python', '-m', 'cbig.CBIG_ukbb_dnn_train', '--phe_dir', phe_dir,
        '--sub_dir', sub_dir, '--data_dir', data_dir, '--icv_dir', icv_dir,
        '--out_dir', out_dir, '--out_subdir', out_subdir, '--inter_dir',
        inter_dir, '--epochs', epochs, '--start_idx', phe_start_idx,
        '--end_idx', phe_end_idx, '--channel_number', *cn, '--dataset',
        dataset, '--batch_size', batch_size, '--lr', lr,
        '--scheduler_decrease', scheduler_decrease, '--output_dim', output_dim,
        '--dropout', dropout
    ])


def test_DNN_outputs(args):
    '''function to test DNN testing function
    Args:
       args: args from command line

    Returns:
        None
     '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    dataset = args.setup + "_test_dataset"

    in_dir = os.path.join(args.base_dir, args.setup, 'data')
    phe_dir = os.path.join(in_dir, args.setup + '_test' + '_final.csv')
    sub_dir = os.path.join(in_dir, args.setup + '_test' + '_subj_list.txt')
    icv_dir = os.path.join(in_dir, args.setup + '_test' + '_icv.csv')
    ukbb_icv_dir = os.path.join(in_dir, args.setup + '_train' + '_icv.csv')
    tra_sub_dir = os.path.join(in_dir,
                               args.setup + '_train' + '_subj_list.txt')
    data_dir = os.path.join(in_dir, 'T1')
    out_dir = os.path.join(res_dir, args.setup + "_test_dataset")
    out_subdir = "basemodel_output"
    model_dir = os.path.join(res_dir,
                             args.setup + "_train_dataset/trained_model_ukbb")
    batch_size = "200"
    output_dim = "3"
    phe_start_idx = "1"
    phe_end_idx = "3"
    c_dim = "2"
    index = "1"

    # generate inputs for mm_stacking
    subprocess.call([
        'python',
        '-m',
        'cbig.CBIG_ukbb_dnn_test',
        '--phe_dir',
        phe_dir,
        '--sub_dir',
        sub_dir,
        '--data_dir',
        data_dir,
        '--out_dir',
        out_dir,
        '--out_subdir',
        out_subdir,
        '--dataset',
        dataset,
        '--batch_size',
        batch_size,
        '--index',
        index,
        '--outputdim',
        output_dim,
        '--across-dataset',
        '--model_dir',
        model_dir,
        '--icv_dir',
        icv_dir,
        '--tra_sub_dir',
        tra_sub_dir,
        '--ukbb_icv_dir',
        ukbb_icv_dir,
        '--start_idx',
        phe_start_idx,
        '--end_idx',
        phe_end_idx,
    ])

    # generate inputs for classical transfer and mm_finetune
    subprocess.call([
        'python',
        '-m'
        'cbig.CBIG_ukbb_dnn_output',
        '--phe_dir',
        phe_dir,
        '--sub_dir',
        sub_dir,
        '--data_dir',
        data_dir,
        '--out_dir',
        out_dir,
        '--out_subdir',
        out_subdir,
        '--dataset',
        dataset,
        '--batch_size',
        batch_size,
        '--index',
        index,
        '--model_dir',
        model_dir,
        '--icv_dir',
        icv_dir,
        '--ukbb_icv_dir',
        ukbb_icv_dir,
        '--tra_sub_dir',
        tra_sub_dir,
        '--c_dim',
        c_dim,
        '--across-dataset',
        '--start_idx',
        phe_start_idx,
        '--end_idx',
        phe_end_idx,
    ])


def test_DNN_classical_transfer(args):
    '''function to test classical transfer learning function
    Args:
        args: args from command line

    Returns:
        None
    '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    dataset = args.setup + "_test_dataset"

    in_dir = os.path.join(args.base_dir, args.setup, 'data')
    phe_dir = os.path.join(in_dir, args.setup + '_test' + '_final.csv')
    sub_dir = os.path.join(in_dir, args.setup + '_test' + '_subj_list.txt')
    icv_dir = os.path.join(in_dir, args.setup + '_test' + '_icv.csv')
    ukbb_icv_dir = os.path.join(in_dir, args.setup + '_train' + '_icv.csv')
    tra_sub_dir = os.path.join(in_dir,
                               args.setup + '_train' + '_subj_list.txt')
    data_dir = os.path.join(res_dir,
                            args.setup + "_test_dataset/basemodel_output")
    out_dir = os.path.join(res_dir, args.setup + "_test_dataset")
    out_subdir = "/classical_transfer"
    inter_dir = os.path.join(out_dir, "output_intermediate")
    model_dir = os.path.join(res_dir,
                             args.setup + "_train_dataset/trained_model_ukbb")
    batch_size = "50"
    index = "1"
    phe_start_idx = "1"
    phe_end_idx = "3"
    epochs = "1"
    ks = ["10", "20", "50", "100", "200"]
    n_rng = "2"
    in_c = "2"
    out_c = "2"

    subprocess.call([
        'python', '-m', 'cbig.CBIG_ukbb_classical_transfer', '--phe_dir',
        phe_dir, '--sub_dir', sub_dir, '--data_dir', data_dir, '--out_dir',
        out_dir, '--out_subdir', out_subdir, '--start_idx', phe_start_idx,
        '--end_idx', phe_end_idx, '--dataset', dataset, '--batch_size',
        batch_size, '--index', index, '--across-dataset', '--model_dir',
        model_dir, '--icv_dir', icv_dir, '--ukbb_icv_dir', ukbb_icv_dir,
        '--tra_sub_dir', tra_sub_dir, '--epochs', epochs, '--ks', *ks,
        '--inter_dir', inter_dir, '--n_rng', n_rng, '--in_c', in_c, '--out_c',
        out_c
    ])


def test_DNN_mm_fintune(args):
    '''function to test meta-matching finetune function
    Args:
        args: args from command line

    Returns:
        None
    '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    dataset = args.setup + "_test_dataset"

    in_dir = os.path.join(args.base_dir, args.setup, 'data')
    phe_dir = os.path.join(in_dir, args.setup + '_test' + '_final.csv')
    sub_dir = os.path.join(in_dir, args.setup + '_test' + '_subj_list.txt')
    icv_dir = os.path.join(in_dir, args.setup + '_test' + '_icv.csv')
    ukbb_icv_dir = os.path.join(in_dir, args.setup + '_train' + '_icv.csv')
    tra_sub_dir = os.path.join(in_dir,
                               args.setup + '_train' + '_subj_list.txt')
    data_dir = os.path.join(res_dir,
                            args.setup + "_test_dataset/basemodel_output")
    out_dir = os.path.join(res_dir, args.setup + "_test_dataset")
    out_subdir = "/mm_finetune"
    inter_dir = os.path.join(out_dir, "output_intermediate")
    model_dir = os.path.join(res_dir,
                             args.setup + "_train_dataset/trained_model_ukbb")
    batch_size = "50"
    index = "1"
    phe_start_idx = "1"
    phe_end_idx = "3"
    epochs = "1"
    ks = ["10", "20", "50", "100", "200"]
    n_rng = "2"
    in_c = "2"
    out_c = "2"

    subprocess.call([
        'python', '-m', 'cbig.CBIG_ukbb_mm_finetune', '--phe_dir', phe_dir,
        '--sub_dir', sub_dir, '--data_dir', data_dir, '--out_dir', out_dir,
        '--out_subdir', out_subdir, '--start_idx', phe_start_idx, '--end_idx',
        phe_end_idx, '--dataset', dataset, '--batch_size', batch_size,
        '--index', index, '--across-dataset', '--model_dir', model_dir,
        '--icv_dir', icv_dir, '--ukbb_icv_dir', ukbb_icv_dir, '--tra_sub_dir',
        tra_sub_dir, '--epochs', epochs, '--ks', *ks, '--inter_dir', inter_dir,
        '--n_rng', n_rng, '--in_c', in_c, '--out_c', out_c
    ])


def test_DNN_mm_stacking(args):
    '''function to test meta-matching stacking function
    Args:
        args: args from command line

    Returns:
        None
    '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    dataset = args.setup + "_test_dataset"

    in_dir = os.path.join(args.base_dir, args.setup, 'data')
    phe_dir = os.path.join(in_dir, args.setup + '_test' + '_final.csv')
    sub_dir = os.path.join(in_dir, args.setup + '_test' + '_subj_list.txt')
    data_dir = os.path.join(res_dir,
                            args.setup + "_test_dataset/basemodel_output")
    out_dir = os.path.join(res_dir, args.setup + "_test_dataset")
    inter_dir = os.path.join(out_dir, "output_intermediate")
    out_subdir = "/mm_stacking"
    phe_start_idx = "1"
    phe_end_idx = "3"
    ks = ["10", "20", "50", "100", "200"]
    n_rng = "2"
    k_limit = "3"
    metric = 'cod'

    subprocess.call([
        'python', '-m', 'cbig.CBIG_ukbb_mm_stacking', '--phe_dir', phe_dir,
        '--sub_dir', sub_dir, '--data_dir', data_dir, '--out_dir', out_dir,
        '--out_subdir', out_subdir, '--start_idx', phe_start_idx, '--end_idx',
        phe_end_idx, '--dataset', dataset, '--across-dataset', '--ks', *ks,
        '--n_rng', n_rng, '--k_limit', k_limit, '--metric', metric,
        '--inter_dir', inter_dir
    ])


def print_results(args):
    '''function to print all results
    Args:
        args: args from command line

    Returns:
        None
    '''

    res_dir = os.path.join(args.base_dir, args.setup, 'results')
    out_dir = os.path.join(res_dir, "examples_test_dataset")

    # print
    print('>>>>>>> Elasticnet <<<<<<<<')
    res = np.load(os.path.join(out_dir,
                               "elasticnet/elasticnet_result_test.npz"),
                  allow_pickle=True)
    print('average correlation:',
          *list(np.mean(np.mean(res.f.meta_cor, axis=2), axis=0)))
    print(
        'average correlation for K = 10       20       50       100      200'
    )
    print('average COD:',
          *list(np.mean(np.mean(res.f.meta_cod, axis=2), axis=0)))
    print('average COD for K = 10        20        50        100       200')

    print('>>>>>>> Directly transfer <<<<<<<<')
    res = np.load(os.path.join(out_dir,
                               "classical_transfer/meta_result_test.npz"),
                  allow_pickle=True)
    print('average correlation:',
          *list(np.mean(np.mean(res.f.meta_cor_tuned, axis=2), axis=0)))
    print(
        'average correlation for K = 10       20       50       100      200'
    )
    print('average COD:',
          *list(np.mean(np.mean(res.f.meta_cod_tuned, axis=2), axis=0)))
    print('average COD for K = 10        20        50        100       200')

    print('>>>>>>> MM finetune <<<<<<<<')
    res = np.load(os.path.join(out_dir, "mm_finetune/meta_result_test.npz"),
                  allow_pickle=True)
    print('average correlation:',
          *list(np.mean(np.mean(res.f.meta_cor_tuned, axis=2), axis=0)))
    print(
        'average correlation for K = 10       20       50       100      200'
    )
    print('average COD:',
          *list(np.mean(np.mean(res.f.meta_cod_tuned, axis=2), axis=0)))
    print('average COD for K = 10        20        50        100       200')

    print('>>>>>>> MM stacking <<<<<<<<')
    res = np.load(os.path.join(out_dir,
                               "mm_stacking/meta_stacking_result_test.npz"),
                  allow_pickle=True)
    print('average correlation:',
          *list(np.mean(np.mean(res.f.meta_cor, axis=2), axis=0)))
    print(
        'average correlation for K = 10       20       50       100      200'
    )
    print('average COD:',
          *list(np.mean(np.mean(res.f.meta_cod, axis=2), axis=0)))
    print('average COD for K = 10        20        50        100       200')


def example_wrapper(args):
    '''example wrapper function
    Args:
        args: args from command line

    Returns:
        None
    '''

    seed = 0
    random.seed(seed)
    np.random.seed(seed)

    # step 1
    generate_data(args)

    # step 2
    test_elasticnet(args)

    # step 3
    test_DNN_MM_training(args)

    # step 4
    test_DNN_outputs(args)

    # step 5
    test_DNN_classical_transfer(args)

    # step 6
    test_DNN_mm_fintune(args)

    # step 7
    test_DNN_mm_stacking(args)

    # step 8
    print_results(args)


if __name__ == "__main__":
    example_wrapper(example_args_parser(setup='examples'))
