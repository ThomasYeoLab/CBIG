"""
Written by Tianchu Zeng and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import argparse
import configparser
import torch
import numpy as np
import os
from CBIG_DELSSOME_pFIC_optimizer import pFIC_Validator, pFIC_Tester, train_help_function
from CBIG_DELSSOME_pFIC_get_EI import get_EI_ratio
from utils.CBIG_DELSSOME_pFIC_utils import convert_csv_to_tensor, convert_mat_to_tensor, set_torch_default


def optimize_body(mode,
                  config_path,
                  seed_nbr,
                  input_dir,
                  output_dir,
                  num_epochs,
                  trials,
                  train_mode,
                  next_epoch=0,
                  seed=None,
                  val_by_dl=False):
    """There are four modes in this function:
    1. train: train pFIC model
    2. validation: validate pFIC model
    3. test: test pFIC model
    4. EI: calculate EI ratio

    Args:
        mode (str): choices from ['train', 'validation', 'test', 'EI']
        config_path (str): config file path
        seed_nbr (int): for distinguish multiple seeds, normally integer from 1 to 5
        input_dir (str): the directory storing input files
        output_dir (str): the directory used to store output files
        num_epochs (int): number of CMA-ES training epochs
        trials (int): Max number of restarting allowed during CMA-ES training when failed
        train_mode (str): choices from ['DELSSOME', 'Euler']
        next_epoch (int, optional): for continual training. Default to 0.
        seed (int, optional): random seed for replication. Default to None and it will be saved out automatically
        val_by_dl (bool, optional): whether to validate by deep learning. Default to False.
            If True, it will use the validation loss from deep learning to select the best parameters.
            If False, it will use the validation loss from Euler method to select the best parameters.
    """

    if not os.path.exists(config_path):
        raise Exception("The config file path does not exist.")

    seed_nbr = int(seed_nbr)
    num_epochs = int(num_epochs)
    trials = int(trials)
    if train_mode == 'DELSSOME':
        dl_epoch_range = np.arange(0, num_epochs)
        euler_epoch_range = []
    elif train_mode == 'Euler':
        dl_epoch_range = []
        euler_epoch_range = np.arange(0, num_epochs)
    else:
        raise Exception("train_mode not valid. Please choose from ['DELSSOME', 'Euler']")
    if seed is not None:
        seed = int(seed)

    config = configparser.ConfigParser()
    config.read(config_path)

    # Read inputs
    myelin = convert_csv_to_tensor(input_dir, 'myelin.csv')
    rsfc_gradient = convert_csv_to_tensor(input_dir, 'rsfc_gradient.csv')
    sc_mat = convert_csv_to_tensor(input_dir, f'SC_{mode}.csv')
    sc_euler = sc_mat / torch.max(sc_mat) * 0.02
    if mode in ['train', 'validation', 'test']:
        fc_emp = convert_csv_to_tensor(input_dir, f'FC_{mode}.csv')
        fcd_cdf_emp = convert_mat_to_tensor(input_dir, f'FCD_CDF_{mode}.mat', 'FCD_CDF')
        fcd_cdf_emp = fcd_cdf_emp / fcd_cdf_emp[-1, 0]

    if mode == 'train':
        save_param_dir = os.path.join(output_dir, f'train/seed{seed_nbr}')
        state = train_help_function(config=config,
                                    myelin=myelin,
                                    RSFC_gradient=rsfc_gradient,
                                    sc_mat=sc_mat,
                                    fc_emp=fc_emp,
                                    fcd_cdf_emp=fcd_cdf_emp,
                                    save_param_dir=save_param_dir,
                                    num_epochs=num_epochs,
                                    dl_epoch_range=dl_epoch_range,
                                    euler_epoch_range=euler_epoch_range,
                                    trials=trials,
                                    next_epoch=next_epoch,
                                    seed=seed)
        return state

    elif mode == 'validation':
        save_param_dir = os.path.join(output_dir, f'train/seed{seed_nbr}')
        val_dir = os.path.join(output_dir, f'validation/seed{seed_nbr}')
        mfm_validator = pFIC_Validator(config, save_param_dir, val_dir, val_by_dl=val_by_dl)
        mfm_validator.val_best_parameters_all_epochs(sc_euler, fc_emp, fcd_cdf_emp, num_epochs, seed=seed)

    elif mode == 'test':
        val_dirs = [os.path.join(output_dir, f'validation/seed{i}') for i in range(seed_nbr, seed_nbr + 1)]
        test_dir = os.path.join(output_dir, f'test/seed{seed_nbr}')
        mfm_tester = pFIC_Tester(config, val_dirs, test_dir, trained_epochs=num_epochs)
        mfm_tester.test(sc_euler, fc_emp, fcd_cdf_emp, seed=seed)

    elif mode == 'EI':
        # select and save the parameter with best validation loss
        val_dirs = [os.path.join(output_dir, f'validation/seed{i}') for i in np.arange(1, seed_nbr + 1)]
        val_best_dir = os.path.join(output_dir, f'val_best/seed{seed_nbr}')
        mfm_tester = pFIC_Tester(config, val_dirs, val_best_dir, trained_epochs=num_epochs)
        mfm_tester.select_best_from_val()

        # Extract the best parameter
        param_path = os.path.join(val_best_dir, 'val_results.pth')
        if not os.path.exists(param_path):
            raise Exception("This group doesn't have validation results.")
        best_parameter = torch.load(param_path)['parameter']

        # Compute and save E/I raio
        EI_save_dir = os.path.join(output_dir, 'EI_ratio')
        EI_save_path = os.path.join(EI_save_dir, f'seed{seed_nbr}.pth')
        os.makedirs(EI_save_dir, exist_ok=True)
        get_EI_ratio(config, save_path=EI_save_path, parameter=best_parameter, sc_euler=sc_euler, seed=seed)

    else:
        raise Exception("Mode not valid.")

    return 0


def optimize_with_args(args, seed_nbr, mode=None):
    if mode is None:
        mode = args.mode
    state = optimize_body(mode=mode,
                          config_path=args.config_path,
                          seed_nbr=seed_nbr,
                          input_dir=args.input_dir,
                          output_dir=args.output_dir,
                          num_epochs=args.num_epochs,
                          trials=args.trials,
                          train_mode=args.train_mode,
                          next_epoch=args.next_epoch,
                          seed=args.seed[seed_nbr - args.start_seed_nbr] if args.seed is not None else None,
                          val_by_dl=args.val_by_dl)
    return state


def optimize_main(args):
    if args.seed is not None:
        assert args.n_seed == len(args.seed)
        # If you specify seed, you need to specify all seeds

    last_seed = args.n_seed + args.start_seed_nbr - 1

    for seed_nbr in range(args.start_seed_nbr, last_seed + 1):

        if args.mode in ['train', 'validation', 'test', 'EI']:
            optimize_with_args(args, seed_nbr)
        elif args.mode == 'train_val_test' or args.mode == 'train_val_EI':
            # Each seed is trained independently; validation selects the best checkpoint;
            # test/EI applies the best checkpoint and saves final metrics.
            state_train = optimize_with_args(args, seed_nbr, mode='train')
            if state_train == 0:
                optimize_with_args(args, seed_nbr, mode='validation')
                if args.mode == 'train_val_test':
                    optimize_with_args(args, seed_nbr, mode='test')
            else:
                print(f"Training failed for seed {seed_nbr}.")

    if args.mode == 'train_val_EI':
        optimize_with_args(args, last_seed, mode='EI')


if __name__ == "__main__":

    parser = argparse.ArgumentParser('pFIC settings')
    parser.add_argument(
        '--mode',
        type=str,
        choices=['train_val_test', 'train_val_EI', 'train', 'validation', 'test', 'EI'],
        default='train',
    )
    parser.add_argument('--config_path',
                        type=str,
                        required=True,
                        help='The path to config file')
    parser.add_argument(
        '--n_seed',
        type=int,
        default=1,
        help=('Number of random initializations (seeds), run one after another in this process. '
              'To use several seeds, prefer one process per seed (--n_seed 1 --start_seed_nbr N) '
              'so they run in parallel on separate CPU cores; see group_level/README.md.'))
    parser.add_argument('--input_dir',
                        type=str,
                        required=True,
                        help='The directory storing input files, including myelin, rsfc_gradient, SC, FC and FCD_CDF')
    parser.add_argument('--output_dir',
                        type=str,
                        required=True,
                        help='The directory storing the output (parameters, losses and E/I ratio)')
    parser.add_argument('--num_epochs', type=int, default=100, help='The total number of CMA-ES training epochs')
    parser.add_argument('--trials',
                        type=int,
                        default=10,
                        help='The number of trials to select new seeds. If seed is specified, force to be 1.')
    parser.add_argument('--train_mode',
                        type=str,
                        choices=['DELSSOME', 'Euler'],
                        default='DELSSOME',
                        help=('The training mode. Choose from DELSSOME (default, the cost predictor) '
                              'or Euler (direct simulation)'))
    parser.add_argument('--next_epoch',
                        type=int,
                        default=0,
                        help='Resume CMA-ES training from the specified epoch index (0-based index)')
    parser.add_argument(
        '--GPU_index',
        type=int,
        default=-1,
        help=('The index of GPU used for training. If no GPU is found, the code will automatically run on CPU.'
              'disabled currently.'))
    parser.add_argument(
        '--seed',
        type=int,
        default=None,
        nargs='+',
        help=('Manually set a random seed for results replication.'
              'If set to None, the function will randomly choose a seed and save it out for future replication.'))
    parser.add_argument('--start_seed_nbr', type=int, default=1, help='The starting seed number')
    parser.add_argument('--val_by_dl',
                        action='store_true',
                        help=('If set, the validation will be done by deep learning. '
                              'If not set, the validation will be done by Euler method.'))

    args = parser.parse_args()

    set_torch_default(-1)
    optimize_main(args)
