#!/usr/bin/env python
# coding: utf-8
'''
### Extract feature importance random forest models for ChenOoi2023_ICCW
This python script contains scripts to extract feature importance from
the random forest models in ChenOoi2023_ICCW.

### Prerequisites
This code assumes that you have:
1. Already finishing training the models and saved them in a .sav format
This code imports the saved model from the .sav format.

Written by Leon Ooi and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

########################################################
# import packages
########################################################
from datetime import datetime
from sklearn.tree import export_text
from sklearn.utils import check_random_state
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
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
# initialize functions for analysis
########################################################
# function to get bootstrap samples from RF tree
def _generate_sample_indices(random_state, n_samples):
    """
    This function regenerates the indices that were sampled
    given a specified random state
    """
    random_instance = check_random_state(random_state)
    sample_indices = random_instance.randint(0, n_samples, n_samples)

    return sample_indices


def _generate_unsampled_indices(random_state, n_samples):
    """
    This function regenerates the indices that were left out
    given a specified random state
    """
    sample_indices = _generate_sample_indices(random_state, n_samples)
    sample_counts = np.bincount(sample_indices, minlength=n_samples)
    unsampled_mask = sample_counts == 0
    indices_range = np.arange(n_samples)
    unsampled_indices = indices_range[unsampled_mask]

    return unsampled_indices


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

# load model
regr = pickle.load(
    open(output_path + '/ABCD_RF_behav_' + str(behav_num) + '.sav', 'rb'))

# save predicted y
y_pred = regr.predict(train_X)
# save to csv
y_pred_csv = pd.DataFrame(np.array(y_pred), columns=['y_pred'])
y_pred_csv.to_csv(
    output_path + '/ABCD_RF_behav_' + str(behav_num) + '_ypred.csv',
    index=False)

########################################################
# start conditional feature importance extraction
########################################################
# setttings for ocndition permutation
perm_amt = 5
corr_thresh = 0.1

# initialize variables
X_sample = train_X.to_numpy(dtype='float32')
n_samples = train_X.shape[0]
n_feats = train_X.shape[1]
n_fivals = np.zeros((n_feats, 1))

# extract permutation score for 3 metrics
fivals_mse = np.zeros((n_feats, 1))
fivals_r2 = np.zeros((n_feats, 1))
fivals_corr = np.zeros((n_feats, 1))

# start permutation
print("---------------------")
print("Start permutation =", datetime.now())
e_count = 0
# iterate over each tree
for estimator in regr.estimators_:
    e_count += 1
    print('\t Estimator:' + str(e_count))
    unsampled_indices = _generate_unsampled_indices(estimator.random_state,
                                                    n_samples)

    # get predictor boundaries
    f_idx = []
    thresh = []
    decision_rules = export_text(estimator)
    print(decision_rules)

    for line in decision_rules.split('\n'):
        for substr in str.split(line):
            if "feature" in substr:
                feature_info = substr.split('_')
                f_idx.append(int(feature_info[1]))
            try:
                thresh.append(float(substr))
            except ValueError:
                pass

    # change into tuple
    fidx_thresh_pairs = list(zip(f_idx, thresh))
    uniq_pairs = list(set(fidx_thresh_pairs))

    # get OOB initial val
    p_estimator = estimator.predict(
        X_sample[unsampled_indices, :], check_input=False)
    p_estimator = p_estimator[:, np.newaxis]
    init_mse = mean_squared_error(train_y[unsampled_indices], p_estimator)
    init_r2 = r2_score(train_y[unsampled_indices], p_estimator)
    init_corr = np.corrcoef(train_y[unsampled_indices],
                            np.squeeze(p_estimator))[0, 1]

    for feat, val in uniq_pairs:  # ignore permutation within site
        feat_idx = feat
        perm_mse = np.zeros((perm_amt, 1))
        perm_r2 = np.zeros((perm_amt, 1))
        perm_corr = np.zeros((perm_amt, 1))

        # find number of groups to permute
        perm_dict = {}
        for idx in range(0, len(uniq_pairs)):
            check_feat_idx = uniq_pairs[idx][0]
            # check if same feature, skip if so
            if feat_idx == check_feat_idx:
                continue
            # check if correlation with other features passes threshold
            check_corr = np.corrcoef(X_sample[:, feat_idx],
                                     X_sample[:, check_feat_idx])[0, 1]
            if abs(check_corr) <= corr_thresh:
                continue
            else:
                # add to dictionary
                if str(check_feat_idx) not in perm_dict:
                    perm_dict[str(check_feat_idx)] = []
                perm_dict[str(check_feat_idx)].append(uniq_pairs[idx][1])

        # check if permutation passes threshold
        for r_seed in range(0, perm_amt):
            Xk_perm = X_sample[unsampled_indices, :]
            g_indices = np.ones((Xk_perm.shape[0], 1))
            g_count = 1
            # randomly permute within each group
            if perm_dict:
                for item in perm_dict:
                    g_count = g_count * (len(perm_dict[item]) + 1)
                    # multiply by number of groups
                    g_indices = g_indices * (len(perm_dict[item]) + 1)
                    sorted_bounds = perm_dict[item]
                    sorted_bounds.sort()
                    for bound in sorted_bounds:
                        # label groups with ascending index based on size
                        g_indices[Xk_perm[:, int(item)] >= bound] -= 1

            # permute within groups
            g_indices = np.squeeze(g_indices)
            for g in range(1, g_count + 1):
                if np.sum(g_indices == g) < 5:
                    print("Warning: Less than 5 samples!")
                bef_perm_idx = np.where(g_indices == g)[0]
                af_perm_idx = np.random.RandomState(
                    seed=r_seed).permutation(bef_perm_idx)
                Xk_perm[bef_perm_idx, feat_idx] = Xk_perm[af_perm_idx,
                                                          feat_idx]

            # find accuracy
            p_estimator = estimator.predict(Xk_perm, check_input=False)
            p_estimator = p_estimator[:, np.newaxis]
            # save permuted feature importance metrics for each tree
            perm_mse[r_seed] = mean_squared_error(train_y[unsampled_indices],
                                                  p_estimator)
            perm_r2[r_seed] = r2_score(train_y[unsampled_indices], p_estimator)
            perm_corr[r_seed] = np.corrcoef(train_y[unsampled_indices],
                                            np.squeeze(p_estimator))[0, 1]
        # compute final feature importance metric
        fivals_mse[feat_idx] += np.mean(perm_mse) - init_mse
        fivals_r2[feat_idx] += np.mean(perm_r2) - init_r2
        fivals_corr[feat_idx] += np.mean(perm_corr) - init_corr
        n_fivals[feat_idx] += 1
# average values
fivals_mse /= n_fivals
fivals_r2 /= n_fivals
fivals_corr /= n_fivals

# set nan to zero
fivals_mse[np.isnan(fivals_mse)] = 0
fivals_r2[np.isnan(fivals_r2)] = 0
fivals_corr[np.isnan(fivals_corr)] = 0

print("End permutation =", datetime.now())
print("---------------------")

########################################################
# extract Haufe feature importance
########################################################
demean_y_pred = y_pred - np.mean(y_pred)
demean_X_train = train_X.to_numpy().T - np.mean(train_X)[:, np.newaxis]

fivals_Haufe = np.dot(demean_X_train, demean_y_pred) / train_X.shape[0]
fivals_Haufe = fivals_Haufe[:, np.newaxis]

########################################################
# save results to csv
########################################################
# save accuracy
acc_csv = pd.DataFrame(
    np.array([
        np.corrcoef(regr.predict(test_X), test_y)[0, 1],
        r2_score(test_y, regr.predict(test_X))
    ])).T
acc_csv.set_axis(['corr', 'COD'], axis=1, inplace=True)
acc_csv.to_csv(
    output_path + '/ABCD_RF_behav_' + str(behav_num) + '_predacc.csv',
    index=False)
# save feature importance
fi_csv = pd.DataFrame(
    np.concatenate((fivals_mse, fivals_r2, fivals_corr, fivals_Haufe), axis=1),
    columns=['fi_mse', 'fi_r2', 'fi_corr', 'fi_Haufe'])
fi_csv.to_csv(
    output_path + '/ABCD_RF_behav_' + str(behav_num) + '_fi.csv', index=False)
