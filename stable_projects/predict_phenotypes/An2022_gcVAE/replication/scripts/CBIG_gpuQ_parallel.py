#!/usr/bin/env python
"""
Written by Lijun An and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""
import os
import argparse
import subprocess


def get_args():
    """
    Setting parameters
    """
    # get user home directory
    home = os.path.expanduser('~')
    storage = os.path.join(home, 'storage')
    parser = argparse.ArgumentParser(description='GPU jobs in parallel')
    parser.add_argument(
        '--cmds',
        type=str,
        required=True,
        help='CMDs to run in parallel on GPU')
    parser.add_argument(
        '--log_names', type=str, required=True, help='Names for log files')
    parser.add_argument(
        '--log_dir',
        type=str,
        default=storage,
        help='Log files dir (default: $home/storage)')
    parser.add_argument(
        '--delimiter',
        type=str,
        default=';',
        help='Splitting cmds & log_names (default: ";")')

    args, _ = parser.parse_known_args()

    return args


def main(args):
    """
    Main function for running jobs in parallel on GPU
    """
    cmd_list = args.cmds.split(args.delimiter)
    lognames_list = args.log_names.split(args.delimiter)
    assert len(cmd_list) == len(lognames_list), \
        'len(cmd_list) ~= len(lognames_list)'
    # use subprocess to run multiple jobs
    processes = []
    for cmd, log_name in zip(cmd_list, lognames_list):
        OUT = open(os.path.join(args.log_dir, log_name + '_out.txt'), 'w')
        ERR = open(os.path.join(args.log_dir, log_name + '_err.txt'), 'w')
        proc = subprocess.Popen(cmd, stdout=OUT, stderr=ERR, shell=True)
        processes.append(proc)
    for proc in processes:
        proc.wait()


if __name__ == '__main__':
    main(get_args())
