# /usr/bin/env python
'''
Written by Kong Xiaolu and CBIG under MIT license:
https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
'''

import os
import numpy as np
import time
import torch
import CBIG_pMFM_basic_functions_main as fc
import warnings


def get_init(myelin_data, gradient_data, highest_order, init_para):
    '''
    This function is implemented to calculate the initial parametrized
    coefficients
    '''

    n_node = myelin_data.shape[0]
    amatrix = np.zeros((n_node, highest_order + 1))
    bmatrix = np.zeros((n_node, highest_order + 1))
    for i in range(highest_order + 1):
        amatrix[:, i] = myelin_data**(i)
        bmatrix[:, i] = gradient_data**(i)
    cmatrix = np.hstack((amatrix, bmatrix[:, 1:highest_order + 1]))
    para = np.linalg.inv(cmatrix.T @ cmatrix) @ cmatrix.T @ init_para
    return para, cmatrix


def CBIG_mfm_SOMA_desikan_main(gpu_index=0, random_seed=1):
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
    Returns:
        None
    '''

    torch.cuda.set_device(gpu_index)
    random_seed_cuda = random_seed
    random_seed_np = random_seed
    torch.manual_seed(random_seed_cuda)
    rng = np.random.Generator(np.random.PCG64(random_seed_np))

    output_path = '../output/step1_training_results/'
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    '''parameter initialization'''
    highest_order = 1
    N = 3 * (2 * highest_order + 1) + 1
    myelin_data = fc.csv_matrix_read('../input/myelin.csv')
    myelin_data = myelin_data[:, 0]
    gradient_data = fc.csv_matrix_read('../input/rsfc_gradient.csv')
    gradient_data = gradient_data[:, 0]
    n_node = myelin_data.shape[0]

    Specimen = np.zeros((68 * 3 + 1, 2))
    Specimen[0:n_node, :] = [0, 1]
    Specimen[n_node:n_node * 2, :] = [0, 0.5]
    Specimen[n_node * 2, :] = [1, 10]
    Specimen[n_node * 2 + 1:, :] = [0.0005, 0.01]
    input_dim = Specimen.shape[0]

    Step = 0.033
    PathLength = 3
    Popsize = 17
    PRT = 0.2
    Migration = 500
    MinDiv = 0.01
    n_dup = 3
    '''create population'''
    Pop_input = rng.standard_normal(size=(input_dim, Popsize)) * np.tile(
        (Specimen[:, 1] - Specimen[:, 0]), [Popsize, 1]).T + np.tile(
            Specimen[:, 0], [Popsize, 1]).T
    Population_w, template_mat = get_init(myelin_data, gradient_data,
                                          highest_order, Pop_input[0:68, :])
    Population_I, template_mat = get_init(myelin_data, gradient_data,
                                          highest_order, Pop_input[68:136, :])
    Population_sigma, template_mat = get_init(
        myelin_data, gradient_data, highest_order, Pop_input[137:205, :])
    Population = np.concatenate(
        (Population_w, Population_I, Pop_input[136:137, :], Population_sigma),
        0)
    arfitness, corr_cost, ks_cost = fc.CBIG_combined_cost_train(
        Pop_input, n_dup)
    '''Migrating'''
    jumpnum = int(PathLength / Step) + 1
    pop_save = np.zeros((N + 3, Migration))

    for gen in range(Migration):

        start_time = time.time()

        leaderindex = np.argmin(arfitness)
        leader = Population[:, leaderindex]
        leader_input = np.zeros(input_dim)
        leader_input[0:68] = template_mat @ leader[0:2 * highest_order + 1]
        leader_input[68:136] = template_mat @ leader[2 * highest_order + 1:2 *
                                                     (2 * highest_order + 1)]
        leader_input[136] = leader[2 * (2 * highest_order + 1)]
        leader_input[137:input_dim] = template_mat @ leader[2 * (
            2 * highest_order + 1) + 1:N]

        np.delete(Population, leaderindex, axis=1)
        costmat = 100 * np.ones((Popsize - 1, jumpnum))
        jumpreal = np.zeros(Popsize - 1)
        jumppop = 0
        popcand = np.zeros((N, jumpnum * (Popsize - 1)))
        pop_input = np.zeros((input_dim, jumpnum * (Popsize - 1)))

        if (np.max(arfitness) -
                np.min(arfitness)) < MinDiv or np.min(arfitness) < 0.4:
            break

        for i in range(Popsize - 1):
            PRTvector = 1 * (rng.standard_normal(N) < PRT)
            realcount = 0
            for t in np.arange(0, PathLength, Step):
                candnew = Population[:, i] + (
                    leader - Population[:, i]) * PRTvector * t
                cand_input = np.zeros(input_dim)
                cand_input[0:68] = template_mat @ candnew[0:2 * highest_order +
                                                          1]
                cand_input[68:136] = template_mat @ candnew[
                    2 * highest_order + 1:2 * (2 * highest_order + 1)]
                cand_input[136] = candnew[2 * (2 * highest_order + 1)]
                cand_input[137:input_dim] = template_mat @ candnew[2 * (
                    2 * highest_order + 1) + 1:N]
                if (cand_input > Specimen[:, 0]).all() & (
                        cand_input < Specimen[:, 1]).all():
                    popcand[:, jumppop] = candnew
                    pop_input[:, jumppop] = cand_input
                    jumppop = jumppop + 1
                    realcount = realcount + 1
            jumpreal[i] = realcount
        popcand = popcand[:, 0:jumppop]
        popcand = np.column_stack((popcand, leader))

        pop_input = pop_input[:, 0:jumppop]
        pop_input = np.column_stack((pop_input, leader_input))

        costM, ks_cost, corr_cost = fc.CBIG_combined_cost_train(
            pop_input, n_dup)

        arfitness[Popsize - 1] = costM[-1]
        cost = costM[0:-1]

        jumpreal = jumpreal.astype(np.int64)
        for i in range(Popsize - 1):
            costmat[i, 0:jumpreal[i]] = cost[np.sum(jumpreal[0:i]):np.
                                             sum(jumpreal[0:i + 1])]
            jumpindex = np.argmin(costmat[i, :])
            Population[:, i] = popcand[:, np.sum(jumpreal[0:i]) + jumpindex]
            arfitness[i] = np.min(costmat[i, :])

        Population[:, Popsize - 1] = leader

        pop_save[0:N, gen] = leader
        pop_save[N, gen] = corr_cost[np.argmin(costM)]
        pop_save[N + 1, gen] = ks_cost[np.argmin(costM)]
        pop_save[N + 2, gen] = np.min(arfitness)

        torch.cuda.empty_cache()

        elapsed = time.time() - start_time
        print('Generation: ', gen + 1)
        print('Cost now: ', np.min(arfitness))
        print('G value: ',
              Population[2 * (2 * highest_order + 1), leaderindex])
        print('Correlation: ', corr_cost[np.argmin(costM)])
        print('KS Cost: ', ks_cost[np.argmin(costM)])
        print('Total number of position: ', corr_cost.shape[0])
        print('Elasped time is: ', elapsed)
        print('****************************************************')

    save_name = [output_path
                 ] + ['random_initialization_',
                      str(random_seed), '.csv']
    np.savetxt(''.join(save_name), pop_save, delimiter=',')


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    for i in range(1, 2):
        CBIG_mfm_SOMA_desikan_main(gpu_index=0, random_seed=i)
