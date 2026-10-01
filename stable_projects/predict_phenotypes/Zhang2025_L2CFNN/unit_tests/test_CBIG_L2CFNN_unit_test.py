#!/usr/bin/env python3
# Written by Chen Zhang and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

import os
import subprocess
import shlex
import shutil
import time
import unittest


def check_job_status(job_name):
    """
    Fetch job status from headnode and return how many jobs with
    "job_name" are still running
    """
    cmd = f'ssh headnode "qstat -f | grep -C 1 `whoami` | grep {job_name} | wc -l"'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def parse_test_log(log_path):
    """
    Parse the unittest log to extract test results and output logs
    """
    test_passed_flag = None
    test_message_lines = []
    capture_message = False

    with open(log_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("Test_passed_flag:"):
                test_passed_flag = line.split(":", 1)[1].strip() == "True"

            elif line.startswith("Test_message:"):
                capture_message = True
                test_message_lines.append(line.split(":", 1)[1].strip())

            elif capture_message:
                # Keep capturing all lines until EOF
                test_message_lines.append(line)

    test_message = "\n".join(test_message_lines)
    return test_passed_flag, test_message


class TestL2C(unittest.TestCase):
    """
    Test data processing, train & evaluation of L2C-FNN,
    L2C-XGBnw, and L2C-XGBw
    """

    def test_L2C(self):

        # --- 1. Setup paths and essential variables ---
        cbig_code_dir = os.environ.get("CBIG_CODE_DIR")
        if not cbig_code_dir:
            print("Error: CBIG_CODE_DIR environment variable not set.")
            return 1

        root_dir = os.path.join(
            cbig_code_dir,
            "stable_projects",
            "predict_phenotypes",
            "Zhang2025_L2CFNN",
        )

        # This path is used for creating log directories locally before submission
        cbig_l2c_root_path = os.path.join(root_dir, "unit_tests")

        print(
            f"INFO: CBIG_L2C_ROOT_PATH for logs (created locally): {cbig_l2c_root_path}"
        )
        os.makedirs(cbig_l2c_root_path, exist_ok=True)

        # --- 2. Define job parameters ---
        log_dir = os.path.join(cbig_l2c_root_path, "job_logs")
        os.makedirs(log_dir, exist_ok=True)

        job_log = os.path.join(log_dir, "unittest.out")
        job_err = os.path.join(log_dir, "unittest.err")
        job_name = "L2CFNN_unittest"

        # --- 3. Construct the command to be executed on the compute node (`inner_cmd`) ---
        # This command includes all necessary environment setup for the job itself.

        # Define PYTHONPATH for the job: $ROOTDIR:$PYTHONPATH
        # Use existing PYTHONPATH from the submission environment as a fallback.
        existing_pythonpath = os.environ.get("PYTHONPATH", "")
        job_pythonpath = (
            f"{root_dir}:{existing_pythonpath}"
            if existing_pythonpath
            else root_dir
        )

        inner_cmd_parts = [
            f"cd {shlex.quote(root_dir)}",
            "source CBIG_init_conda",
            "conda activate CBIG_Zhang2025_py38",  # Ensure conda is initialized in the job's shell
            # Use ${HOME} for the job's execution environment
            "export PYTHONPYCACHEPREFIX=${HOME}/.cache/Python",  # Escaped for the final ssh command string
            f"export PYTHONPATH={shlex.quote(job_pythonpath)}",
            f"export CBIG_L2C_ROOT_PATH={shlex.quote(cbig_l2c_root_path)}",  # Exporting for the job's context
            "python -m unit_tests.L2CFNN_unit_test",  # Your unit test command
        ]
        inner_cmd = "; ".join(inner_cmd_parts)

        submit_cmd = (
            f"ssh headnode \"echo '{inner_cmd}' | "
            "qsub -V -q gpuQ-L40S "
            f"-l walltime=00:10:00 "
            f"-l select=1:ncpus=1:mem=10G:ngpus=1 "
            f"-m ae "
            f"-N {job_name} "
            f'-e {job_err} -o {job_log}"'
        )

        # --- 4. Execute the submission command ---
        try:
            # Using shell=True because ssh_command_str is a complex command string
            # intended for shell interpretation on the local machine (to run ssh).
            result = subprocess.run(
                submit_cmd, shell=True, capture_output=True, text=True
            )
            print("INFO: Job submission command executed successfully.")
            if result.stdout:
                print("INFO: Submission script stdout:", result.stdout)
            if (
                result.stderr
            ):  # Should be empty if successful, but good to check
                print("INFO: Submission script stderr:", result.stderr)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: Job submission failed with exit code {e.returncode}"
            )
            if e.stdout:
                print("ERROR: Submission script stdout:", e.stdout)
            if e.stderr:
                print("ERROR: Submission script stderr:", e.stderr)
        except FileNotFoundError:
            print(
                "ERROR: `ssh` command not found. Please ensure it's in your PATH."
            )

        print(
            ">>>>>> Unittest submitted to scheduler, waiting for job to finish. <<<<<<"
        )
        while check_job_status(job_name):
            time.sleep(10)

        # Parse job log to get test status and per-test results
        test_passed_flag, test_message = parse_test_log(job_log)

        # Remove job logs (intermediate files)
        if os.path.exists(log_dir):
            shutil.rmtree(log_dir)
        else:
            print(f"Folder not found: {log_dir}")

        self.assertTrue(test_passed_flag, test_message)
        print(test_message)


if __name__ == "__main__":
    unittest.main()
