#!/bin/sh
# Written by Chen Zhang and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

import argparse
import random
import shutil
import subprocess

import numpy as np

from config import global_config


class ExampleConfigs:
    """
    Example configuration class with default values.
    Users can override values by passing arguments to the constructor.
    Note:
        All values are stored as strings since they are passed directly to
        subprocess.call() as command-line arguments, which require string inputs.
    """

    def __init__(
        self,
        nb_folds="20",
        site="ADNI",
        fold="0",
        start_fold="0",
        seed="3",
        epoch="1",
        trial="10",
        batch_size="512",
        input_dim="101",
        nb_measures="8",
        num_folds="1",
    ):
        self.nb_folds = nb_folds
        self.site = site
        self.fold = fold
        self.start_fold = start_fold
        self.seed = seed
        self.epoch = epoch
        self.trial = trial
        self.batch_size = batch_size
        self.input_dim = input_dim
        self.nb_measures = nb_measures
        self.num_folds = num_folds


def get_default_config():
    return ExampleConfigs()


example_configs = get_default_config()


def test_data_generation():
    """
    function to generate dummy data
    """
    subprocess.call(
        [
            "python",
            "-m",
            "examples.generate_dummy_data",
            "--site",
            example_configs.site,
            "--seed",
            example_configs.seed,
        ]
    )
    print(">>>>>>>>> Dummy data generation finished")


def test_data_processing():
    """
    function to process test data
    """

    # Step 1: split generation
    subprocess.call(
        [
            "python",
            "-m",
            "data_processing.step1_fold_split",
            "--nb_folds",
            example_configs.nb_folds,
            "--site",
            example_configs.site,
            "--need_split",
        ]
    )
    print(">>>>>>>>> Data processing step 1 finished")

    # Step 2: baseline generation
    subprocess.call(
        [
            "python",
            "-m",
            "data_processing.step2_baseline_gen",
            "--site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> Data processing step 2 finished")

    # Step 3: data per fold generation
    subprocess.call(
        [
            "python",
            "-m",
            "data_processing.step3_fold_gen",
            "--site",
            example_configs.site,
            "--need_split",
        ]
    )
    print(">>>>>>>>> Data processing step 3 finished")

    # Step 4: L2C feature generation (only 1 fold for demonstration purpose)
    subprocess.call(
        [
            "python",
            "-m",
            "data_processing.step4_frog_gen",
            "--site",
            example_configs.site,
            "--fold",
            example_configs.fold,
            "--need_split",
        ]
    )
    print(">>>>>>>>> Data processing step 4 finished")

    # Step 5: process L2C features for L2C-FNN (only 1 fold)
    subprocess.call(
        [
            "python",
            "-m",
            "data_processing.step5_gauss_data_process",
            "--site",
            example_configs.site,
            "--fold",
            example_configs.fold,
            "--train_site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> Data processing step 5 finished")

    # Step 6: filter out unused ground truth data
    subprocess.call(
        [
            "python",
            "-m",
            "data_processing.step6_gt_filter",
            "--site",
            example_configs.site,
            "--need_split",
        ]
    )
    print(">>>>>>>>> Data processing step 6 finished")

    # Step 7: MinimalRNN data process (only 1 fold)
    p1 = subprocess.Popen(
        [
            "python",
            "-m",
            "data_processing.step7_mrnn_gen",
            "--site",
            example_configs.site,
            "--strategy",
            "model",
            "--fold",
            example_configs.fold,
            "--in_domain",
            "--validation",
        ]
    )

    p2 = subprocess.Popen(
        [
            "python",
            "-m",
            "data_processing.step7_mrnn_gen",
            "--site",
            example_configs.site,
            "--strategy",
            "model",
            "--fold",
            example_configs.fold,
            "--in_domain",
        ]
    )

    # Wait for all Step 7 jobs to finish
    p1.wait()
    p2.wait()
    print(">>>>>>>>> Data processing step 8 finished")


def test_L2C_FNN():
    """
    function to test the performance of L2C_FNN (in-domain)
    """

    # Step 1: L2C-FNN training with hyperparameter tuning (only 1 epoch, 10 trials)
    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_FNN.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--epochs",
            example_configs.epoch,
            "--trials",
            example_configs.trial,
            "--batch_size",
            example_configs.batch_size,
            "--seed",
            example_configs.seed,
            "--input_dim",
            example_configs.input_dim,
            "--site",
            example_configs.site,
        ]
    )

    # Step 2: remove checkpoints and logs that are not best hyperparameter setup
    subprocess.call(
        [
            "python",
            "-m",
            "models.model_utils",
            "--util",
            "clear_temp",
            "--model",
            "L2C_FNN",
            "--site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> L2C-FNN training finished")

    # Step 3: in_domain evaluation
    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_FNN.evaluation",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--batch_size",
            example_configs.batch_size,
            "--input_dim",
            example_configs.input_dim,
            "--in_domain",
            "--train_site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> L2C-FNN evaluation finished")


def test_L2C_XGBw():
    """
    function to test the performance of L2C_XGBw (in-domain)
    """

    # Step 1: L2C-XGBw training with hyperparameter tuning (10 boost_rounds, 10 trials)
    p1 = subprocess.Popen(
        [
            "python",
            "-m",
            "models.L2C_XGBw.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--debug",
            "--target",
            "mmse",
        ]
    )

    p2 = subprocess.Popen(
        [
            "python",
            "-m",
            "models.L2C_XGBw.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--debug",
            "--target",
            "vent",
        ]
    )

    p3 = subprocess.Popen(
        [
            "python",
            "-m",
            "models.L2C_XGBw.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--debug",
            "--target",
            "clin",
        ]
    )

    # Wait for all 3 training processes to finish
    p1.wait()
    p2.wait()
    p3.wait()

    # Step 2: remove checkpoints and logs that are not best hyperparameter setup
    subprocess.call(
        [
            "python",
            "-m",
            "models.model_utils",
            "--util",
            "clear_temp",
            "--model",
            "L2C_XGBw",
            "--site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> L2C-XGBw training finished")

    # Step 3: in_domain evaluation
    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_XGBw.eval",
            "--fold",
            example_configs.fold,
            "--site",
            example_configs.site,
            "--in_domain",
            "--train_site",
            example_configs.site,
            "--target",
            "mmse",
        ]
    )

    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_XGBw.eval",
            "--fold",
            example_configs.fold,
            "--site",
            example_configs.site,
            "--in_domain",
            "--train_site",
            example_configs.site,
            "--target",
            "clin",
        ]
    )

    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_XGBw.eval",
            "--fold",
            example_configs.fold,
            "--site",
            example_configs.site,
            "--in_domain",
            "--train_site",
            example_configs.site,
            "--target",
            "vent",
        ]
    )
    print(">>>>>>>>> L2C-XGBw evaluation finished")


def test_L2C_XGBnw():
    """
    function to test the performance of L2C_XGBnw (in-domain)
    """

    # Step 1: L2C-XGBnw training with hyperparameter tuning (10 boost_rounds, 10 trials)
    p1 = subprocess.Popen(
        [
            "python",
            "-m",
            "models.L2C_XGBnw.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--debug",
            "--target",
            "mmse",
        ]
    )

    p2 = subprocess.Popen(
        [
            "python",
            "-m",
            "models.L2C_XGBnw.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--debug",
            "--target",
            "vent",
        ]
    )

    p3 = subprocess.Popen(
        [
            "python",
            "-m",
            "models.L2C_XGBnw.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--debug",
            "--target",
            "clin",
        ]
    )

    # Wait for all 3 training processes to finish
    p1.wait()
    p2.wait()
    p3.wait()

    # Step 2: remove checkpoints and logs that are not best hyperparameter setup
    subprocess.call(
        [
            "python",
            "-m",
            "models.model_utils",
            "--util",
            "clear_temp",
            "--model",
            "L2C_XGBnw",
            "--site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> L2C-XGBnw training finished")

    # Step 3: in_domain evaluation
    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_XGBnw.eval",
            "--fold",
            example_configs.fold,
            "--site",
            example_configs.site,
            "--in_domain",
            "--train_site",
            example_configs.site,
            "--target",
            "mmse",
        ]
    )

    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_XGBnw.eval",
            "--fold",
            example_configs.fold,
            "--site",
            example_configs.site,
            "--in_domain",
            "--train_site",
            example_configs.site,
            "--target",
            "clin",
        ]
    )

    subprocess.call(
        [
            "python",
            "-m",
            "models.L2C_XGBnw.eval",
            "--fold",
            example_configs.fold,
            "--site",
            example_configs.site,
            "--in_domain",
            "--train_site",
            example_configs.site,
            "--target",
            "vent",
        ]
    )
    print(">>>>>>>>> L2C-XGBnw evaluation finished")


def test_MinimalRNN():
    """
    function to test the performance of MinimalRNN (in-domain)
    """

    # Step 1: MinimalRNN training with hyperparameter tuning (only 1 epoch, 10 trials)
    subprocess.call(
        [
            "python",
            "-m",
            "models.MinimalRNN.train",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--epochs",
            example_configs.epoch,
            "--trials",
            example_configs.trial,
            "--batch_size",
            example_configs.batch_size,
            "--seed",
            example_configs.seed,
            "--nb_measures",
            example_configs.nb_measures,
            "--site",
            example_configs.site,
        ]
    )

    # Step 2: remove checkpoints and logs that are not best hyperparameter setup
    subprocess.call(
        [
            "python",
            "-m",
            "models.model_utils",
            "--util",
            "clear_temp",
            "--model",
            "MinimalRNN",
            "--site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> MinimalRNN training finished")

    # Step 3: in_domain evaluation
    subprocess.call(
        [
            "python",
            "-m",
            "models.MinimalRNN.evaluation",
            "--fold",
            example_configs.fold,
            "--start_fold",
            example_configs.start_fold,
            "--site",
            example_configs.site,
            "--batch_size",
            example_configs.batch_size,
            "--in_domain",
            "--train_site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> MinimalRNN evaluation finished")


def test_AD_Map():
    """
    function to test the performance of AD_Map (in-domain)
    """

    # Step 1: AD_Map training with hyperparameter tuning (n_iter_base=10, 10 trials)
    subprocess.call(
        [
            "python",
            "-m",
            "models.AD_Map.train",
            "--fold",
            example_configs.fold,
            "--trials",
            example_configs.trial,
            "--seed",
            example_configs.seed,
            "--site",
            example_configs.site,
            "--debug",
        ]
    )

    # Step 2: remove checkpoints and logs that are not best hyperparameter setup
    subprocess.call(
        [
            "python",
            "-m",
            "models.model_utils",
            "--util",
            "clear_temp",
            "--model",
            "AD_Map",
            "--site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> AD_Map training finished")

    # Step 3: in_domain evaluation
    subprocess.call(
        [
            "python",
            "-m",
            "models.AD_Map.evaluation",
            "--fold",
            example_configs.fold,
            "--in_domain",
            "--site",
            example_configs.site,
            "--train_site",
            example_configs.site,
        ]
    )
    print(">>>>>>>>> AD_Map evaluation finished")


def check_results():
    """
    function to compare example run prediction against reference results
    """

    subprocess.call(
        [
            "python",
            "-m",
            "utils.check_reference_results",
            "--models",
            "L2C_FNN",
            "L2C_XGBw",
            "L2C_XGBnw",
            "AD_Map",
            "MinimalRNN",
            "--site",
            example_configs.site,
            "--num_folds",
            example_configs.num_folds,
        ]
    )
    print(">>>>>>>>> Results comparison finished")


def clean_up(folder_list=None):
    """
    Deletes specified folders located in the same parent directory
    as this script. If no folder list is provided, uses a default list.

    Args:
        folder_list (list, optional): List of folder names to delete.
    """

    # Get the directory where the current script is located
    # script_dir = Path(__file__).parent.resolve() # Use resolve() for absolute path
    script_dir = global_config.root_path

    # List of folder names to delete relative to the script's directory
    default_folders = ["checkpoints", "data", "predictions", "raw_data"]
    folders_to_delete = folder_list if folder_list is not None else default_folders

    print(f"Script directory: {script_dir}")
    print("Attempting to delete the following folders:")
    print(folders_to_delete)

    for folder_name in folders_to_delete:
        # Construct the full path to the target folder
        folder_path = script_dir / folder_name

        # Check if the path exists and is actually a directory
        if folder_path.is_dir():
            print(f"Deleting folder: {folder_path} ... ", end="")
            try:
                # Recursively delete the directory and all its contents
                shutil.rmtree(folder_path)
                print("Done.")
            except OSError as e:
                # Handle potential errors (e.g., permission denied)
                print(f"Error deleting {folder_path}: {e.strerror}")
            except Exception as e:
                # Catch other potential exceptions
                print(
                    f"An unexpected error occurred while deleting {folder_path}: {e}"
                )
        elif folder_path.exists():
            # Path exists but is not a directory (it's a file)
            print(
                f"Skipping '{folder_name}': It exists but is a file, not a folder."
            )
        else:
            # Path does not exist
            print(f"Skipping '{folder_name}': Folder not found.")

    print("Cleanup process finished.")


def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--step", type=int, required=True)
    return parser


def example_wrapper(args):
    """
    example wrapper function
    """

    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    if args.step == 1:
        # step 0
        print(">>>>>> Start data generation <<<<<<")
        test_data_generation()

        # step 1
        print(">>>>>> Start data processing <<<<<<")
        test_data_processing()

        # step 2
        print(">>>>>> Start L2C-FNN train and evaluation <<<<<<")
        test_L2C_FNN()

        # step 3
        print(">>>>>> Start L2C-XGBw train and evaluation <<<<<<")
        test_L2C_XGBw()

        # step 4
        print(">>>>>> Start L2C-XGBnw train and evaluation <<<<<<")
        test_L2C_XGBnw()

        # step 5
        print(">>>>>> Start MinimalRNN train and evaluation <<<<<<")
        test_MinimalRNN()
    elif args.step == 2:
        # step 6
        print(">>>>>> Start AD_Map train and evaluation <<<<<<")
        test_AD_Map()

        # step 7
        print(">>>>>> Compare results against reference <<<<<<")
        check_results()

        # step 8
        print(">>>>>> Cleanup generated files <<<<<<")
        clean_up()
    else:
        print(f"Example step {args.step} is not defined, please choose 1 or 2")


if __name__ == "__main__":
    args = get_arg_parser().parse_args()
    example_wrapper(args)
