#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Naren Wulan and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import unittest
import numpy as np
import random
from utils.check_reference_results import check_reference_results
from examples.CBIG_MMT1_example import example_args_parser, generate_data, \
 test_elasticnet, test_DNN_MM_training, test_DNN_outputs, \
 test_DNN_classical_transfer, test_DNN_mm_fintune, test_DNN_mm_stacking


class TestMM(unittest.TestCase):
    """
    Test the training of elasticnet, classical transfer learning,
    meta-matching finetune and stacking
    """

    def test_training(self):
        args = example_args_parser(setup='unit_tests')

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
        check_reference_results(args)


if __name__ == '__main__':
    unittest.main()
