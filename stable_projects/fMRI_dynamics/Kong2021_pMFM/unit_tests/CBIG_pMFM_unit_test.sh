#!/bin/bash
# this function is the unit test
# Written by Kong Xiaolu and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

cd ../examples/scripts
source activate pMFM
python CBIG_pMFM_example_check_results.py
