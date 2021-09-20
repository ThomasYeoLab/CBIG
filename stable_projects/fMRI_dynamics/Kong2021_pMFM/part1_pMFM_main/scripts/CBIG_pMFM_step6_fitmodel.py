# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import scipy.io as sio
from sklearn.mixture import GaussianMixture


def CBIG_pMFM_fitmodel_empirical():
    '''
    This function is to implement the FCD distribution fit
    for empirical data
    Returns:        None
    '''

    FCD_data_raw = sio.loadmat(
        '../output/step5_STDFCD_results/STD_FCD_empirical_rundata.mat')
    FCD_mean_all = FCD_data_raw['FCD_emp_allrun']
    file_num = FCD_mean_all.shape[0]

    output_path = '../output/step6_SWSTD_state/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    bic_all = np.zeros((file_num, 3))
    mid_p = np.zeros((file_num, 1))

    for i in range(file_num):
        FCD_mean = np.expand_dims(FCD_mean_all[i, :], 1)

        gm1 = GaussianMixture(n_components=1, n_init=1)
        gm1.fit(FCD_mean)

        gm2 = GaussianMixture(n_components=2, n_init=10)
        gm2.fit(FCD_mean)

        bic1 = gm1.bic(FCD_mean)
        bic2 = gm2.bic(FCD_mean)

        bic_all[i, 0] = bic1
        bic_all[i, 1] = bic2
        bic_all[i, 2] = np.argmin(bic_all[i, 0:2])

        if bic_all[i, 2] == 1:
            mu1 = gm2.means_[0]
            sigma1 = gm2.covariances_[0, 0, 0]
            mu2 = gm2.means_[1]
            sigma2 = gm2.covariances_[1, 0, 0]
            delta = gm2.weights_[0]

            a = 1 / sigma1**2 - 1 / sigma2**2
            b = 2 * (mu2 / sigma2**2 - mu1 / sigma1**2)
            c = mu1**2 / sigma1**2 - mu2**2 / sigma2**2 - 2 * np.log(
                delta / (1 - delta) * sigma2 / sigma1)

            mid_p1 = (-b + np.sqrt(b**2 - 4 * a * c)) / 2 / a
            mid_p2 = (-b - np.sqrt(b**2 - 4 * a * c)) / 2 / a
            if mu1 < mid_p1 < mu2 or mu1 > mid_p1 > mu2:
                mid_p[i, 0] = mid_p1
            else:
                mid_p[i, 0] = mid_p2

        print(i)

    threshold = dict()
    threshold['th_all_emp'] = mid_p
    threshold['bic_all_emp'] = bic_all
    sio.savemat(output_path + 'FCD_threshold_empirical.mat', threshold)


def CBIG_pMFM_fitmodel_simulated():
    '''
        This function is to implement the FCD distribution fit
        for simulated data
        Returns:        None
        '''

    FCD_data_raw = sio.loadmat(
        '../output/step5_STDFCD_results/STD_FCD_simulated_rundata.mat')
    FCD_mean_all = FCD_data_raw['FCD_sim_allrun']
    file_num = FCD_mean_all.shape[0]

    output_path = '../output/step6_SWSTD_state/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    bic_all = np.zeros((file_num, 3))
    mid_p = np.zeros((file_num, 1))

    for i in range(file_num):
        FCD_mean = np.expand_dims(FCD_mean_all[i, :], 1)

        gm1 = GaussianMixture(n_components=1, n_init=1)
        gm1.fit(FCD_mean)

        gm2 = GaussianMixture(n_components=2, n_init=10)
        gm2.fit(FCD_mean)

        bic1 = gm1.bic(FCD_mean)
        bic2 = gm2.bic(FCD_mean)

        bic_all[i, 0] = bic1
        bic_all[i, 1] = bic2
        bic_all[i, 2] = np.argmin(bic_all[i, 0:2])

        if bic_all[i, 2] == 1:
            mu1 = gm2.means_[0]
            sigma1 = gm2.covariances_[0, 0, 0]
            mu2 = gm2.means_[1]
            sigma2 = gm2.covariances_[1, 0, 0]
            delta = gm2.weights_[0]

            a = 1 / sigma1**2 - 1 / sigma2**2
            b = 2 * (mu2 / sigma2**2 - mu1 / sigma1**2)
            c = mu1**2 / sigma1**2 - mu2**2 / sigma2**2 - 2 * np.log(
                delta / (1 - delta) * sigma2 / sigma1)

            mid_p1 = (-b + np.sqrt(b**2 - 4 * a * c)) / 2 / a
            mid_p2 = (-b - np.sqrt(b**2 - 4 * a * c)) / 2 / a
            if mu1 < mid_p1 < mu2 or mu1 > mid_p1 > mu2:
                mid_p[i, 0] = mid_p1
            else:
                mid_p[i, 0] = mid_p2

        print(i)

    threshold = dict()
    threshold['th_all_sim'] = mid_p
    threshold['bic_all_sim'] = bic_all
    sio.savemat(output_path + 'FCD_threshold_simulated.mat', threshold)


if __name__ == "__main__":
    CBIG_pMFM_fitmodel_empirical()
    CBIG_pMFM_fitmodel_simulated()
