# /usr/bin/env python
"""
Written by Shaoshi Zhang and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
"""

import time
import os
import shutil
import sys
from multiprocessing import Pool
from multiprocessing import set_start_method
import configparser
import numpy as np
import torch
import CBIG_pFIC_misc as misc


def CBIG_mfm_optimization_main(parser, process_id, random_seed):
    """

    This function optimizes model parameters of the neural mass model by
    running CMAES to minimize the cost function, which consists of a FC
    correlation cost, a FC L1 cost and a FCD KS cost.
    This approach optimizes the 9 linear coefficients and a global scaling
    constant (in total, 10 parameters). This function utilizes python
    multiprocessing to maximize the GPU usage.
    -Input:
        -parser: a parsed configuration file specifying model parameters
            (see example.ini)
        -process_id: a number used as an identifier for a process
        -random_seed: a number used to set a seed for RNG
    -Output:
        -training_iteration.txt: a log recording the current CMAES epoch
        -cost.txt: a log recoding the training progress
        -training.csv: optimized parameters and associated costs. Each row
            corresponds to 1 CMAES epoch, the first 10 numbers are the
            linear coeff of wEE, linear coeff of wEI, G and linear coeff of
            sigma. The last 4 numbers of FC correlation loss, FC L1 loss,
            FCD KS loss and total loss.
    [Terminology] 1.'Model parameter/regional parameter/parameter' refers to
                    regional wEE, wEI or sigma.
                  2.'Linear coefficients' refer to the linear coefficients of
                    model parameter. For example wEE = a + b*myelin +
                    c*RSFC_gradient, here, a, b, and c are the linear
                    coefficients. Together with G (global scaling constant),
                    these are the variables that are directly optimized during
                    CMAES.

    """

    output_path = parser['training']['output_path']
    output_file = parser['training']['output_file'] + '_' + str(
        process_id) + '.csv'
    FCD = parser['training']['FCD_training']
    SC = parser['training']['SC_training']
    FC = parser['training']['FC_training']
    rE_min = float(parser['system']['rE_min'])
    rE_max = float(parser['system']['rE_max'])
    FC_corr_weight = float(parser['training']['FC_corr_weight'])
    FC_L1_weight = float(parser['training']['FC_L1_weight'])
    FCD_weight = float(parser['training']['FCD_KS_weight'])
    myelin = parser['training']['myelin_training']
    rsfc_gradient = parser['training']['RSFC_gradient_training']
    gpu_number = int(parser['system']['GPU_index'])
    if gpu_number != -1:
        torch.cuda.set_device(gpu_number)
    gpu_name = torch.cuda.get_device_name(0)

    # Set random seed for numpy and pytorch
    random_seed_cuda = random_seed
    random_seed_np = random_seed
    torch.manual_seed(random_seed_cuda)
    rng = np.random.Generator(np.random.PCG64(random_seed_np))

    # Initialize input parameters
    myelin_data = misc.csv_matrix_read(myelin)
    gradient_data = misc.csv_matrix_read(rsfc_gradient)

    num_myelin_component = myelin_data.shape[1]
    num_gradient_component = gradient_data.shape[1]
    # number of linear coefficients and G
    N = 3 * (num_myelin_component + num_gradient_component + 1) + 1
    # number of linear coeffs associated with each parameter
    N_p = num_myelin_component + num_gradient_component + 1
    n_node = myelin_data.shape[0]
    dim = n_node * 3 + 1  # total number of parameters

    search_range = misc.set_search_range(parser, misc.csv_matrix_read(FC))
    init_para = rng.uniform(0, 1, dim) * (
        search_range[:, 1] - search_range[:, 0]) + search_range[:, 0]
    start_point_wEE, concat_matrix = \
        get_init(myelin_data, gradient_data, init_para[0:n_node])
    start_point_wEI, concat_matrix = \
        get_init(myelin_data, gradient_data, init_para[n_node:n_node * 2])
    start_point_sigma, concat_matrix = \
        get_init(myelin_data, gradient_data, init_para[n_node * 2 + 1:dim])

    # Initialize CMAES children processes
    xmean = np.zeros(N)
    xmean[0:N_p] = start_point_wEE
    xmean[1] = start_point_wEE[1] / 2
    xmean[N_p:2 * N_p] = start_point_wEI
    xmean[N_p + 1] = start_point_wEI[1] / 2
    xmean[2 * N_p] = init_para[2 * n_node]
    xmean[2 * N_p + 1:N] = start_point_sigma

    # Initialize optimization hyper-parameters
    cov = 0.25  # coordinate wise standard deviation (step size)
    # total number of CMAES epochs
    num_epoch = int(parser['training']['num_epoch'])
    n_dup = 5

    # CMAES parameters setting
    Lambda = 100  # number of children processes
    mu = 10
    weights = np.log(mu + 1 / 2) - np.log(np.arange(1, mu + 1))
    weights = weights / np.sum(weights)
    mueff = np.sum(weights)**2 / np.sum(weights**2)

    # Strategy parameter setting: adaptation
    cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N)
    cs = (mueff + 2) / (N + mueff + 5)
    c1 = 2 / ((N + 1.3)**2 + mueff)
    cmu = np.minimum(1 - c1,
                     2 * (mueff - 2 + 1 / mueff) / ((N + 2)**2 + mueff))
    damps = 1 + 2 * np.maximum(0, np.sqrt((mueff - 1) / (N + 1)) - 1) + cs

    # Initializing dynamic strategy parameters and constants
    pc = np.zeros(N)
    ps = np.zeros(N)
    B = np.eye(N)
    D = np.ones(N)
    D[0:N_p] = start_point_wEE[0] / 2
    D[N_p:2 * N_p] = start_point_wEI[0] / 2
    D[2 * N_p] = 0.4
    D[2 * N_p + 1:N] = 0.001 / 2

    C = np.dot(np.dot(B, np.diag(np.power(D, 2))), B.T)
    invsqrtC = np.dot(np.dot(B, np.diag(np.power(D, -1))), B.T)
    chiN = N**0.5 * (1 - 1 / (4 * N) + 1 / (21 * N ^ 2))

    # Evolution loop
    epoch = 0
    arx = np.zeros([N, Lambda])
    input_para = np.zeros((dim, Lambda))
    xmin = np.zeros(N + 4)

    # spatial constraint refers to the assumption that spatial distribution
    # of wEI and myelin should have a negative correlation, while spatial
    # distribution of wEI and rsfc gradient should have a positive correlation
    spatial_constraint = int(parser['training']['spatial_constraint'])
    with open(output_path + output_file, 'a') as f:
        while epoch < num_epoch:
            iteration_log = open(
                output_path + 'training_iteration_' + str(random_seed) +
                '.txt', 'w')
            iteration_log.write(str(epoch) + '\n')
            iteration_log.write(gpu_name)
            iteration_log.close()

            start_time = time.time()

            # Generate children processes
            arx[:, 0] = xmean
            j = 0
            infinite_loop_count = 0
            while j < Lambda:
                arx[:, j] = xmean + cov * np.dot(B,
                                                 (D * rng.standard_normal(N)))
                input_para[0:n_node, j] = concat_matrix @ arx[0:N_p, j]
                input_para[n_node:2 *
                           n_node, j] = concat_matrix @ arx[N_p:2 * N_p, j]
                wei_myelin_corr = np.corrcoef(
                    np.squeeze(myelin_data, 1),
                    input_para[n_node:2 * n_node, j])
                wei_myelin_corr = wei_myelin_corr[0, 1]
                wei_rsfc_corr = np.corrcoef(
                    np.squeeze(gradient_data, 1),
                    input_para[n_node:2 * n_node, j])
                wei_rsfc_corr = wei_rsfc_corr[0, 1]
                input_para[2 * n_node:2 * n_node + 1, j] = arx[2 * N_p, j]
                input_para[2 * n_node +
                           1:dim, j] = concat_matrix @ arx[2 * N_p + 1:N, j]

                if spatial_constraint:
                    if (input_para[:, j] < search_range[:, 0]).any() or (
                            input_para[:, j] > search_range[:, 1]).any() \
                            or wei_rsfc_corr < 0 or wei_myelin_corr > 0:
                        j = j - 1
                        infinite_loop_count += 1
                        if infinite_loop_count > 200000:
                            iteration_log = open(
                                output_path + 'training_iteration_' +
                                str(random_seed) + '.txt', 'w')
                            iteration_log.write(str(epoch) + ': Infinite Loop')
                            iteration_log.close()
                            return
                    j = j + 1
                else:
                    if (input_para[:, j] < search_range[:, 0]).any() or (
                            input_para[:, j] > search_range[:, 1]).any():
                        j = j - 1
                        infinite_loop_count += 1
                        if infinite_loop_count > 200000:
                            iteration_log = open(
                                output_path + 'training_iteration_' +
                                str(random_seed) + '.txt', 'w')
                            iteration_log.write(str(epoch) + ': Infinite Loop')
                            iteration_log.close()
                            return
                    j = j + 1

            print("[Process ID: " + str(process_id) + "] CMA-ES epoch " +
                  str(epoch))
            total_cost, fc_corr_cost, fc_L1_cost, fcd_cost, bold, r_E = \
                misc.CBIG_combined_cost_train(input_para, n_dup, SC, FC,
                                              FCD, rE_min, rE_max,
                                              [FC_corr_weight,
                                               FC_L1_weight], FCD_weight,
                                              parser)
            epoch = epoch + 1

            # Sort by total cost and compute weighted mean
            arfitsort = np.sort(total_cost)
            arindex = np.argsort(total_cost)
            xold = xmean
            xmean = np.dot(arx[:, arindex[0:mu]], weights)
            xshow = xmean - xold

            # Cumulation
            ps = (1 - cs) * ps + np.sqrt(cs * (2 - cs) * mueff) * np.dot(
                invsqrtC, xshow) / cov
            hsig = (np.linalg.norm(ps) / np.sqrt(1 - pow(
                (1 - cs), (2 * epoch))) / chiN < (1.4 + 2 / (N + 1))) * 1
            pc = (1 - cc) * pc + hsig * np.sqrt(cc *
                                                (2 - cc) * mueff) * xshow / cov

            # Adapting covariance matrix C
            artmp = (1 / cov) * (
                arx[:, arindex[0:mu]] - np.tile(xold, [mu, 1]).T)
            C = (1 - c1 - cmu) * C + c1 * (
                        np.outer(pc, pc) + (1 - hsig) * cc * (2 - cc) * C) + \
                cmu * np.dot(artmp, np.dot(np.diag(weights), artmp.T))

            # Adapting step size
            cov = cov * np.exp((cs / damps) * (np.linalg.norm(ps) / chiN - 1))

            # Decomposition
            if 1 > 1 / (c1 + cmu) / N / 10:
                C = np.triu(C, k=1) + np.triu(C).T
                D, B = np.linalg.eigh(C)
                D = D.real
                B = B.real
                D = np.sqrt(D)
                invsqrtC = np.dot(B, np.dot(np.diag(D**(-1)), B.T))

            # Monitoring the evolution status
            cost_log = open(output_path + 'cost_' + str(random_seed) + '.txt',
                            'w')
            print('******** Generation: ' + str(epoch) + ' ********')
            cost_log.write('******** Generation: ' + str(epoch) + ' ********' +
                           '\n')
            print('Mean total cost: ', np.mean(arfitsort[0:mu]))

            xmin[0:N] = arx[:, arindex[0]]
            xmin[N] = fc_corr_cost[arindex[0]]
            xmin[N + 1] = fc_L1_cost[arindex[0]]
            xmin[N + 2] = fcd_cost[arindex[0]]
            xmin[N + 3] = np.min(total_cost)
            xmin_save = np.reshape(xmin, (-1, N + 4))
            np.savetxt(f, xmin_save, delimiter=',')
            print('Cov: ', cov)
            cost_log.write('Cov: ' + str(cov) + '\n')
            print('Best parameter set: ', arindex[0])
            cost_log.write('Best parameter set: ' + str(arindex[0]) + '\n')
            print('Best total cost: ', np.min(total_cost))
            cost_log.write('Best total cost: ' + str(np.min(total_cost)) +
                           '\n')
            print('FC correlation cost: ', fc_corr_cost[arindex[0]])
            cost_log.write('FC correlation cost: ' +
                           str(fc_corr_cost[arindex[0]]) + '\n')
            print('FC L1 cost: ', fc_L1_cost[arindex[0]])
            cost_log.write('FC L1 cost: ' + str(fc_L1_cost[arindex[0]]) + '\n')
            print('FCD KS statistics cost: ', fcd_cost[arindex[0]])
            cost_log.write('FCD KS statistics cost: ' +
                           str(fcd_cost[arindex[0]]) + '\n')
            cost_log.write('rE min max: ' + str(rE_min) + ' ,' + str(rE_max) +
                           '\n')
            cost_log.write('Cost function weights: ' + str(FC_corr_weight) +
                           ' ,' + str(FC_L1_weight) + ' ,' + str(FCD_weight) +
                           '\n')

            elapsed_time = time.time() - start_time
            print('Time taken: ', elapsed_time)
            cost_log.write('Time taken: ' + str(elapsed_time) + '\n')
            print('******************************************')
            cost_log.write('******************************************')
            cost_log.close()


def get_init(myelin, gradient, init_para):
    """

    This function apply GLM to an initial set of regional parameters to get
    the linear coefficients.
    -Input:
        -myelin: #ROI x 1 T1w/T2w ratio
        -gradient: #ROI x 1 RSFC gradient
        -init_para: #ROI x 1 regional model parameters (wEE, wEI or sigma)
    -Output:
        -linear_coeff: the linear coefficients used for optimization
        -concat_matrix: a horizontally concatenated matrix consisted of a bias
            term, myelin and RSFC gradient data

    """

    n_node = myelin.shape[0]
    # bias term, myelin, RSFC gradient
    concat_matrix = np.vstack((np.ones(n_node), myelin.T, gradient.T)).T
    linear_coeff = np.linalg.inv(
        concat_matrix.T @ concat_matrix) @ concat_matrix.T @ init_para
    return linear_coeff, concat_matrix


if __name__ == "__main__":
    config = sys.argv[1]
    parser = configparser.ConfigParser()
    parser.read(config)
    num_thread = int(parser['system']['num_thread'])
    random_seed = eval(parser['system']['random_seed'])
    if len(random_seed) != num_thread:
        raise Exception("The number of multiprocessing threads does not "
                        "match the number of random seeds!")

    output_path = parser['training']['output_path']
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    else:
        shutil.rmtree(output_path)
        os.makedirs(output_path)

    input_args = ()
    for i in range(num_thread):
        random_seed_thread = random_seed[i]
        input_args_thread = [parser, i + 1, random_seed_thread]
        input_args = input_args + (input_args_thread, )

    set_start_method("spawn")
    p = Pool(num_thread)
    p.starmap(CBIG_mfm_optimization_main, input_args)
