# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import time
import torch
import CBIG_pMFM_basic_functions as fc
import warnings


def CBIG_mfm_optimization_desikan_main(gpu_index=0, random_seed=1):
    '''
    This function is to implement the optimization processes of mean
    field model.
    The objective function is the summation of FC correlation cost and
    FCD KS statistics cost.
    The optimization process is highly automatic and generate 500
    candidate parameter sets for
    main results.

    Args:
        gpu_index:      index of gpu used for optimization
        random_seed:    random seed for optimization
        output_path:    output directory for saving optimized model
        parameters
    Returns:
        None
    '''

    # Setting random seed and GPU
    torch.cuda.set_device(gpu_index)
    random_seed_cuda = random_seed
    random_seed_np = random_seed
    torch.manual_seed(random_seed_cuda)
    rng = np.random.Generator(np.random.PCG64(random_seed_np))

    # Create output folders
    output_path = '../output/step1_training_results/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # Initializing input parameters
    n_node = 68
    N = n_node * 3 + 1

    search_range = np.zeros((N, 2))
    search_range[0:n_node, :] = [0, 1]
    search_range[n_node:2 * n_node, :] = [0, 0.5]
    search_range[2 * n_node, :] = [1, 10]
    search_range[2 * n_node + 1:, :] = [0.0005, 0.01]
    init_para = rng.uniform(0, 1, N) * (
        search_range[:, 1] - search_range[:, 0]) + search_range[:, 0]
    xmean = init_para.copy()

    # Initializing optimization hyper-parameters
    sigma = 0.15
    sigmaS = 0.15
    stoppoint = 0.3
    maxloop = 500
    n_dup = 3

    # CMA-ES parameters setting
    Lambda = 500
    mu = 40
    weights = np.log(mu + 1 / 2) - np.log(np.arange(1, mu + 1))
    weights = weights / np.sum(weights)
    mueff = 1 / np.sum(weights**2)

    # Strategy parameter setting: adaptation
    cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N)
    cs = (mueff + 2) / (N + mueff + 5)
    c1 = 2 / ((N + 1.3)**2 + mueff)
    cmu = np.minimum(1 - c1,
                     2 * (mueff - 2 + 1 / mueff) / ((N + 2)**2 + mueff))
    damps = 1 + 2 * np.maximum(0, np.sqrt((mueff - 1) / (N + 1)) - 1) + cs

    # Initializing dynamic strategy parameters and constants'''
    pc = np.zeros(N)
    ps = np.zeros(N)
    B = np.eye(N)
    D = np.zeros(N)
    D[0:n_node] = search_range[0, 1] / 2
    D[n_node:2 * n_node] = search_range[n_node, 1] / 2
    D[2 * n_node] = 0.4
    D[2 * n_node + 1:3 * n_node + 1] = search_range[n_node * 2 + 1, 1] / 2
    C = np.dot(np.dot(B, np.diag(np.power(D, 2))), B.T)
    invsqrtC = np.dot(np.dot(B, np.diag(np.power(D, -1))), B.T)
    chiN = N**0.5 * (1 - 1 / (4 * N) + 1 / (21 * N ^ 2))

    # Evolution loop
    countloop = 0
    arx = np.zeros([N, Lambda])
    xmin = np.zeros([N + 3, maxloop])
    stop_count = 0
    while countloop < maxloop:

        start_time = time.time()

        # Generating lambda offspring
        arx[:, 0] = xmean
        j = 0
        while j < Lambda:
            arx[:, j] = xmean + sigma * np.dot(B, (D * rng.standard_normal(N)))
            if (arx[:, j] < search_range[:, 0]).any() or (
                    arx[:, j] > search_range[:, 1]).any():
                arx[:, j] = arx[:, j] + 2 * (
                    arx[:, j] <
                    search_range[:, 0]) * (search_range[:, 0] - arx[:, j])
                arx[:, j] = arx[:, j] + 2 * (
                    arx[:, j] >
                    search_range[:, 1]) * (search_range[:, 1] - arx[:, j])
            j = j + 1

        # Calculating costs of offspring
        total_cost, fc_cost, fcd_cost = fc.CBIG_combined_cost_train(arx, n_dup)
        countloop = countloop + 1

        # Sort by total cost and compute weighted mean
        arfitsort = np.sort(total_cost)
        arindex = np.argsort(total_cost)
        xold = xmean
        xmean = np.dot(arx[:, arindex[0:mu]], weights)
        xshow = xmean - xold

        # Cumulation
        ps = (1 - cs) * ps + np.sqrt(cs * (2 - cs) * mueff) * np.dot(
            invsqrtC, xshow) / sigma
        hsig = (np.linalg.norm(ps) / np.sqrt(1 - (1 - cs)**
                                             (2 * countloop)) / chiN <
                (1.4 + 2 / (N + 1))) * 1
        pc = (1 - cc) * pc + hsig * np.sqrt(cc *
                                            (2 - cc) * mueff) * xshow / sigma

        # Adapting covariance matrix C
        artmp = (1 / sigma) * (
            arx[:, arindex[0:mu]] - np.tile(xold, [mu, 1]).T)
        C = (1 - c1 - cmu) * C + c1 * (
            np.outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C) + cmu * np.dot(
                artmp, np.dot(np.diag(weights), artmp.T))

        # Adapting step size
        sigma = sigma * np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1))
        sigma = min(sigma, sigmaS)

        # Decomposition
        if 1 > 1 / (c1 + cmu) / N / 10:
            C = np.triu(C, k=1) + np.triu(C).T
            if np.isnan(C).any():
                break
            D, B = np.linalg.eigh(C)
            D = D.real
            B = B.real
            D = np.sqrt(D)
            invsqrtC = np.dot(B, np.dot(np.diag(D**(-1)), B.T))

        # Monitoring the evolution status
        ps_norm = np.linalg.norm(ps)
        print('******** Generation: ' + str(countloop) + ' ********')
        print('Norm of P-sigma: ', ps_norm)
        print('The mean of total cost: ', np.mean(arfitsort[0:mu]))
        print('Sigma: ', sigma)

        xmin[0:N, countloop - 1] = arx[:, arindex[0]]
        xmin[N, countloop - 1] = fc_cost[arindex[0]]
        xmin[N + 1, countloop - 1] = fcd_cost[arindex[0]]
        xmin[N + 2, countloop - 1] = np.min(total_cost)
        print('Best total cost: ', np.min(total_cost))
        print('FC correlation cost: ', fc_cost[arindex[0]])
        print('FCD KS statistics cost: ', fcd_cost[arindex[0]])

        elapsed_time = time.time() - start_time
        print('Elapsed time for this evolution is : ', elapsed_time)
        print('******************************************')

        # break
        if arfitsort[0] < stoppoint and ps_norm < 11:
            stop_count = stop_count + 1
        if stop_count >= 5 or sigma < 0.001:
            break

    save_name = [output_path] + ['random_seed_', str(random_seed), '.csv']
    np.savetxt(''.join(save_name), xmin, delimiter=',')


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=UserWarning)
    for i in range(1, 11):
        CBIG_mfm_optimization_desikan_main()
