import os
import time
import numba
import click
import pickle
import numpy as np
from typing import Tuple, List, Any
from random import random, normalvariate, choices

import utils


@numba.jit
def calculate_score(matrix: np.ndarray,
                    events: np.ndarray,
                    target_result: np.ndarray,
                    target_err: np.ndarray,
                    covariance_matrix: np.ndarray,
                    data_without_mec: np.ndarray,
                    normalize: bool=False
                    ) -> Tuple[float, float]:
    cross_sections = np.matmul(matrix, events)
    nuwro_result = cross_sections + data_without_mec

    if normalize:
        nuwro_result *= np.sum(target_result) / np.sum(nuwro_result)

    chi2_cov = utils.chi2(target_result, nuwro_result,
                          covariance_matrix)

    chi2_std = utils.chi2_standard(target_result, nuwro_result,
                                   target_err)

    return [chi2_cov, chi2_std]


def mutate(matrix: np.ndarray,
           mutation_probability: float,
           smoothness: float
           ) -> np.ndarray:
    # iterate over allowed elements
    matrix = matrix.copy().reshape(24, 24)

    for i in range(0, 24):
        for j in range(i, 24):
            # decide if element will be changed
            if random() < mutation_probability:
                # choose mutation type
                if random() < 0.5:
                    matrix[i][j] = np.clip(matrix[i][j] * normalvariate(1, 1),
                                           a_min=0.1, a_max=None)
                else:
                    matrix[i][j] = np.clip(matrix[i][j] + normalvariate(0, 0.5),
                                           a_min=0.1, a_max=None)
                if smoothness > 0:
                    smooth_low = 1 - smoothness
                    smooth_high = 1 / smooth_low
                    neighbours = [matrix[(i + 1) % 24][j],
                                  matrix[i - 1][j],
                                  matrix[i][(j + 1) % 24],
                                  matrix[i][j - 1]]

                    new_neighbours = [x for x in neighbours if x > 0]

                    max_v = smooth_high * min(new_neighbours)
                    min_v = smooth_low * max(new_neighbours)
                    matrix[i][j] = np.clip(matrix[i][j],
                                           a_min=min_v, a_max=max_v)

    return matrix.reshape(24 * 24)


def cross(matrices: List[List[Any]],
          mutation_probability: float,
          smoothness: float
          ) -> np.ndarray:
    child = matrices[0][1].copy().reshape((24, 24))
    parent_b = matrices[1][1].reshape((24, 24))

    probability_o = matrices[0][0] / (matrices[0][0] + matrices[1][0])

    if smoothness > 0:
        smooth_low = 1 - smoothness
        smooth_high = 1 / smooth_low

    for i in range(0, 24):
        for j in range(i, 24):
            if random() > probability_o:
                child[i][j] = parent_b[i][j]
                if smoothness > 0:
                    neighbours = [child[(i + 1) % 24][j],
                                  child[i - 1][j],
                                  child[i][(j + 1) % 24],
                                  child[i][j - 1]]

                    new_neighbours = [x for x in neighbours if x > 0]

                    max_v = smooth_high * min(new_neighbours)
                    min_v = smooth_low * max(new_neighbours)

                    child[i][j] = np.clip(child[i][j],
                                          a_min=min_v, a_max=max_v)

    child.shape = (24 * 24)

    return mutate(child, mutation_probability, smoothness)


def new_population(population_with_score: List[Tuple[Any, np.ndarray]],
                   mutation_probability: float,
                   smoothness: float
                   ) -> List[np.ndarray]:
    chi2, old_population = [], []

    for c, p in population_with_score:
        chi2.append(np.sum(c[0]))
        old_population.append([0, p])

    reduction = max(chi2)

    for i, v in enumerate(chi2):
        chi2[i] = abs(v - reduction)
        old_population[i][0] = chi2[i]

    new_population = [population_with_score[i][1] for i in range(5)]
    for _ in range(5, len(population_with_score)):
        new_population.append(cross(choices(old_population,
                                            weights=chi2, k=2),
                                    mutation_probability, smoothness))

    return new_population


@click.command()
@click.option('--iterations', default=20, type=int)
@click.option('--plot', is_flag=True)
@click.option('--show_matrix', is_flag=True)
@click.option('--data_path')
@click.option('--experiment', required=True)
@click.option('--population_size', default=10, type=int)
@click.option('--mutation_probability', default=0.05, type=float)
@click.option('--smoothness', default=0, type=float)
@click.option('--initial_matrix', default="ones", type=str)
@click.option('--eval_type', default="chi2_cov", type=str)
@click.option('--breed', is_flag=True)
@click.option('--breed_disable', default=0, type=int)
@click.option('--save', is_flag=True)
@click.option('--shape_covariance', is_flag=True)
@click.option('--normalize', is_flag=True)
@click.option('--pseudoinverse', is_flag=True)
def genetic(iterations: int,
            plot: bool,
            show_matrix: bool,
            data_path: str,
            experiment: str,
            population_size: int,
            mutation_probability: float,
            initial_matrix: str,
            smoothness: float,
            eval_type: str,
            breed: bool,
            breed_disable: int,
            save: bool,
            shape_covariance: bool,
            normalize: bool,
            pseudoinverse: bool
            ) -> None:

    if initial_matrix == 'ones':
        original_matrix = np.ones((24, 24))
    else:
        with open(initial_matrix, 'rb') as f:
            original_matrix = pickle.load(f).reshape((24, 24))

    for i in range(24):
        for j in range(i + 1, 24):
            original_matrix[j][i] = 0
    original_matrix = original_matrix.reshape(24 * 24)

    experiment_info = utils.load_experimental_data(experiment)
    target_res, target_err = experiment_info[:2]
    covariance_matrix, data_without_mec = experiment_info[2:4]
    trans_left, long_left = experiment_info[4:6]

    if not data_path:
        data_path = os.path.join("Data", experiment)
    data = utils.import_data_from_path(data_path)

    target_shape = target_res.shape

    events = data.events.reshape((24 * 24, -1))
    target_res = target_res.reshape(-1)
    target_err = target_err.reshape(-1)
    data_without_mec = data_without_mec.reshape(-1)

    if "minerva" in experiment or experiment == "optimize_test":
        events = utils.events_to_cross_section(events,
                                               data.nof_events,
                                               data.cross_section,
                                               target_shape,
                                               trans_left,
                                               long_left)
    if shape_covariance:
        if experiment == 'minerva_neutrino':
            import covariance_matrix
            covariance = covariance_matrix.covariance
            # _, covariance = covariance_matrix.inv(covariance)
        elif experiment == 'minerva_antyneutrino':
            from minerva_antyneutrino import minerva_antyneutrino_covariance
            covariance = minerva_antyneutrino_covariance
        shape_matrix, _, _ = utils.shape_matrix_extraction(target_res,
                                                           covariance)
        # print(np.linalg.det(shape_matrix))
        if pseudoinverse:
            covariance_matrix = np.linalg.pinv(shape_matrix)
        else:
            covariance_matrix = np.linalg.inv(shape_matrix)

    initial_chi2 = calculate_score(original_matrix,
                                   events,
                                   target_res,
                                   target_err,
                                   covariance_matrix,
                                   data_without_mec,
                                   normalize)
    print(f"Chi2 for chosen initial_matrix: {initial_chi2}")

    population = [original_matrix]
    population += [mutate(original_matrix, mutation_probability, smoothness)
                   for _ in range(1, population_size)]

    if eval_type == 'chi2_cov':
        key = lambda m: m[0][0]
    elif eval_type == 'chi2_std':
        key = lambda m: m[0][1]
    elif eval_type == 'chi2_both':
        key = lambda m: m[0][0] + m[0][1]

    smoothness = smoothness if smoothness < 1 and smoothness > 0 else 0

    start = time.time()
    try:
        for i_no in range(iterations):
            chi2 = [calculate_score(population[i],
                                    events,
                                    target_res,
                                    target_err,
                                    covariance_matrix,
                                    data_without_mec,
                                    normalize)
                    for i in range(population_size)]

            population_with_score = list(zip(chi2, population))

            population_with_score.sort(key=key,
                                       reverse=False)

            print(i_no, population_with_score[0][0])

            if breed_disable:
                if i_no > breed_disable:
                    breed = False

            if breed:
                population = new_population(population_with_score,
                                            mutation_probability,
                                            smoothness)
            else:
                population[0] = population_with_score[0][1]
                for i in range(1, population_size):
                    population[i] = mutate(population_with_score[0][1],
                                           mutation_probability,
                                           smoothness)
    except KeyboardInterrupt:
        print("Manually stopped computation")

    elapsed = time.time() - start
    minutes = elapsed / 60
    seconds = int(elapsed) % 60
    print(f"Calculated in: {int(minutes)}m{seconds}s")

    if save:
        save_chi2 = population_with_score[0][0]
        config = {
            'experiment': experiment,
            'iterations': i_no,
            'data_path': data_path,
            'population_size': population_size,
            'mutation_probability': mutation_probability,
            'initial_matrix': initial_matrix,
            'smoothness': smoothness if smoothness else 1,
            'eval_type': eval_type,
            'breed': breed,
            'chi2_cov': save_chi2[0],
            'chi2_std': save_chi2[1]
        }
        name = f'scale_cov{save_chi2[0]:.2f}_std{save_chi2[1]:.2f}.pkl'
        with open(os.path.join('Matrices', name), 'wb') as file:
            pickle.dump(population_with_score[0][1].reshape(24, 24), file)
            pickle.dump(config, file)
            print(f"Scaling matrix saved to {name}")

    if plot or show_matrix:
        from matplotlib import pyplot

    if plot:
        scaled_data = np.dot(population_with_score[0][1],
                             events).reshape(12, 13)
        scaled_data = utils.events_to_cross_section(scaled_data,
                                                    data.nof_events,
                                                    data.cross_section,
                                                    target_res.shape)

        for i in range(12):
            pyplot.hlines(scaled_data[i] + utils.cross_sections_without_mec[i],
                          xmin=utils.trans_nu_left[:-1],
                          xmax=utils.trans_nu_left[1:], colors='g')
            pyplot.hlines(utils.target_result[i],
                          xmin=utils.trans_nu_left[:-1],
                          xmax=utils.trans_nu_left[1:])
            pyplot.legend(['Scaled NuWro', "MINERvA"])
            pyplot.show()

    if show_matrix:
        utils.show_matrix(population_with_score[0][1].reshape(24, 24))

if __name__ == '__main__':
    genetic()
