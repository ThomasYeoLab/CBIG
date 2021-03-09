# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import csv
import os
import numpy as np
import CBIG_pMFM_parameter_estimation_example as example
import warnings


def csv_matrix_read(filename):
    '''
    This function is used to read csv file into a numpy array
    Args:
        filename:  input csv file
    Returns:
        out_array: output numpy array
    '''

    csv_file = open(filename, "r")
    read_handle = csv.reader(csv_file)
    out_list = []
    R = 0
    for row in read_handle:
        out_list.append([])
        for col in row:
            out_list[R].append(float(col))
        R = R + 1
    out_array = np.array(out_list)
    csv_file.close()
    return out_array


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    print('Start toy example execution')
    example.CBIG_mfm_optimization_desikan_main()
    print('Finish toy example execution')

    print('Start results comparison')
    expected_results = csv_matrix_read('../ref_output/expected_output.csv')
    example_results = csv_matrix_read('../output/example_output.csv')
    os.remove('../output/example_output.csv')
    os.rmdir('../output/')
    if np.linalg.norm(example_results - expected_results) > 0.0001:
        raise ValueError('Unit test failed.')
    else:
        print('Unit test passed.')
