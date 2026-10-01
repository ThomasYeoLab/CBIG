#!/bin/sh
# Written by Chen Zhang and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

import argparse
import random
import subprocess

import numpy as np
from examples.CBIG_L2CFNN_example import (
    clean_up,
    get_default_config,
    test_data_generation,
    test_data_processing,
    test_L2C_FNN,
    test_L2C_XGBnw,
    test_L2C_XGBw,
)

example_configs = get_default_config()


def check_results():
    """
    function to compare example run prediction against reference results
    """

    result = subprocess.run(
        [
            "python",
            "-m",
            "utils.check_reference_results",
            "--models",
            "L2C_FNN",
            "L2C_XGBw",
            "L2C_XGBnw",
            "--site",
            example_configs.site,
            "--num_folds",
            example_configs.num_folds,
        ],
        capture_output=True,
        text=True,
    )

    # Print full output for logging
    # print(result.stdout)
    print("Error output: ", result.stderr)

    # Raise an exception if the return code indicates failure
    return result.returncode == 0, result.stdout


def get_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    return parser


def main():
    """
    Main function calling each sub-module of the test
    and returning the compare_reference results
    """
    args = get_arg_parser().parse_args()

    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)

    # step 0
    test_data_generation()

    # step 1
    test_data_processing()

    # step 2
    test_L2C_FNN()

    # step 3
    test_L2C_XGBw()

    # step 4
    test_L2C_XGBnw()

    # step 5
    print(">>>>>> Compare results against reference <<<<<<")
    passed, msg = check_results()

    # step 6
    print(">>>>>> Cleanup generated files <<<<<<")
    clean_up(['output'])

    print(f"Test_passed_flag: {passed}")
    print(f"Test_message: {msg}")


if __name__ == "__main__":
    main()
