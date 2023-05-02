#!/usr/bin/env python
# coding: utf-8
'''
### Train random forest models for ChenOoi2023_ICCW
This python script contains scripts to train the random forest
models in ChenOoi2023_ICCW.

### Prerequisites
This code assumes that you have:
1. Saved the FC matrix in csv format
2. Already run KRR
This code imports the FC matrix from a csv format, and uses the
splits generated during KRR.

### Inputs
fold_num: An integer of the fold number to run regression for
behav_num: An integer of the behavior number to run regression for
curr_sample: An integer of the sample size to run regression for
results_dir = A string of the output path

### Outputs
sav_file: A sav file of the saved random forest model

Written by Leon Ooi and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

########################################################
# import packages
########################################################
from datetime import datetime
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import scipy.io as sio
import pickle
import sys
import os

########################################################
# initialize paths and inputs
########################################################
# required inputs
fold_num = int(sys.argv[1])
behav_num = int(sys.argv[2])
curr_sample = int(sys.argv[3])
results_dir = sys.argv[4]

# define output and input directories: MODIFY HERE IF NEEDED
krr_dir = results_dir + '/KRR/'
output_dir = results_dir + '/RF/'
fc_csv = os.getenv('CBIG_REPDATA_DIR') + \
 '/stable_projects/predict_phenotypes' + '/ChenOoi2023_ICCW/input' \
 + '/ICCW_5260_FC.csv'

# custom paths to required folders in krr_dir
y_regressed = krr_dir + str(curr_sample) + '/y/fold_' + str(
    fold_num) + '/y_regress_all_score.mat'
folds = krr_dir + str(curr_sample) + '/no_relative_5_fold_sub_list.mat'

# create output folder
output_folder = 'behav_' + str(behav_num) + '/rng0/' + 'fold_' + str(fold_num)
output_path = os.path.join(output_dir, str(curr_sample), output_folder)
if not os.path.exists(output_path):
    os.makedirs(output_path)

# print inputs for reference in log file
print("---------------------")
print("fold_num:" + str(fold_num))
print("behav_num:" + str(behav_num))
print("curr_sample:" + str(curr_sample))
print("Generating results in:" + output_dir)
print("---------------------")

########################################################
# start processing
########################################################
# load files
print("Load files...")
fold_y = sio.loadmat(y_regressed)
all_fold = sio.loadmat(folds)

# get variables from files
behav_num_tmp = behav_num
curr_y = fold_y['y_resid'][:, (behav_num_tmp - 1)]

# pad array to full size
missing_vals = 5260 - len(curr_y)
curr_y = np.pad(curr_y, (0, missing_vals))
fc = pd.read_csv(fc_csv)

# get fold subjects
train_idx = [bool(i) for i in all_fold['sub_fold'][fold_num - 1][0][1] == 0]
test_idx = [bool(i) for i in all_fold['sub_fold'][fold_num - 1][0][1] == 1]

# normalize fc
fc_norm = (fc - np.mean(fc, axis=0)) / np.std(fc, axis=0)

# split into train and test
test_y = curr_y[test_idx]
train_y = curr_y[train_idx]
test_X = fc_norm.iloc[test_idx, :]
train_X = fc_norm.iloc[train_idx, :]

########################################################
# start training
########################################################
# settings
num_trees = 100
c_depth = 4

print("---------------------")
print("Settings: Trees =", num_trees, ", Depth =", c_depth)
print("Start time =", datetime.now())
# train
regr = RandomForestRegressor(
    n_estimators=num_trees, max_depth=c_depth, random_state=0, oob_score=True)
regr.fit(train_X, train_y)
print("Train Acc (corr):", np.corrcoef(regr.predict(train_X), train_y)[0, 1])
print("Test Acc (corr):", np.corrcoef(regr.predict(test_X), test_y)[0, 1])
print("Test Acc (COD):", r2_score(test_y, regr.predict(test_X)))
print("End time =", datetime.now())
print("---------------------")

########################################################
# save model
########################################################
pickle.dump(
    regr, open(output_path + '/ABCD_RF_behav_' + str(behav_num) + '.sav',
               'wb'))
