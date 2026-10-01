#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
"""
Written by Shaoshi Zhang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
import os
import time
import shutil
import unittest
import subprocess
import numpy as np
import CBIG_pFIC_misc as misc

CBIG_CODE_DIR = os.environ.get('CBIG_CODE_DIR')
BASE_DIR = os.path.join(CBIG_CODE_DIR, 'stable_projects', 'fMRI_dynamics',
                        'Zhang2024_pFIC')
UNIT_TESTS_DIR = os.path.join(BASE_DIR, 'unit_tests')
OUTPUT_DIR = os.path.join(UNIT_TESTS_DIR, 'output')
REFERENCE_OUTPUT_DIR = os.path.join(UNIT_TESTS_DIR, 'ref_results')
JOB_NAME = 'pFIC_unittest'


def check_job_status():
    cmd = ('ssh headnode " qstat -f | grep -C 1 `whoami` | grep '
           f'{JOB_NAME} | wc -l"')
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())


class CBIG_pFIC_unit_test(unittest.TestCase):
    def tearDown(self):
        # Remove generated output only if the test passed, so failed runs
        # keep their output around for debugging.
        test_passed = all(
            exc_info is None for _, exc_info in self._outcome.errors)
        if test_passed:
            shutil.rmtree(OUTPUT_DIR, ignore_errors=True)

    def test_pFIC(self):
        config = os.path.join(UNIT_TESTS_DIR, 'CBIG_pFIC_unit_test_config.ini')
        training_script = os.path.join(UNIT_TESTS_DIR,
                                       'CBIG_pFIC_unit_test_training.py')
        transpose_script = os.path.join(
            UNIT_TESTS_DIR, 'CBIG_pFIC_unit_test_transpose_matrix.py')

        cmd = (
            f"cd {BASE_DIR}; "
            "source CBIG_init_conda; conda activate CBIG_pFIC; "
            f"python {training_script} {config}; "
            f"python {transpose_script} {config}; "
            "conda deactivate")
        submit_cmd = (
            'ssh headnode " $CBIG_CODE_DIR/setup/CBIG_pbsubmit '
            f'-cmd \'{cmd}\' -walltime 1:00:00 -mem 5G -ncpus 1 -ngpus 1 '
            f'-name {JOB_NAME}"')
        result = subprocess.run(
            submit_cmd, shell=True, capture_output=True, text=True)
        print('INFO: submission stdout:', result.stdout)
        if result.stderr:
            print('INFO: submission stderr:', result.stderr)

        while check_job_status():
            time.sleep(60)

        user_result = misc.csv_matrix_read(
            os.path.join(OUTPUT_DIR, 'training', 'training_1.csv'))
        ref_result = misc.csv_matrix_read(
            os.path.join(REFERENCE_OUTPUT_DIR, 'training', 'training_1.csv'))
        diff = np.sum(np.abs(user_result - ref_result))
        self.assertLess(
            diff, 1e-6,
            'pFIC training results differ too much from the reference '
            f'outputs (sum abs diff = {diff}).')


if __name__ == '__main__':
    unittest.main()
