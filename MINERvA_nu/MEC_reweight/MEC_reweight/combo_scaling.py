import os
import time
import numba
import click
import pickle
import numpy as np
from typing import Tuple, List, Any

import utils
from genetic import calculate_score, mutate, cross, new_population


class Experiment(object):

    def __init__(self, name, monte_carlo, experimental, bins):
        self.name = name
        self.monte_carlo = monte_carlo.reshape((24 * 24, -1))
        self.result = experimental[0].reshape(-1)
        self.error = experimental[1].reshape(-1)
        self.covariance = experimental[2]
        self.data_without_mec = experimental[3].reshape(-1)
        self.bins = bins
        if 'T2K' in name:
            self.monte_carlo *= 12 / 13
            self.data_without_mec *= 12 / 13


class Combined_score(object):

    def __init__(self, experiments, norm_penalty=0.0):
        self.experiments = experiments
        self.norm_penalty = norm_penalty

    def __call__(self, matrix: np.ndarray,
                 normalize_by_bins: bool,
                 normalize_by_cross_section: bool):
        score = []
        for experiment in self.experiments:
            if normalize_by_cross_section and 'minerva' in experiment.name:
                monte_carlo = np.sum(matrix @ experiment.monte_carlo)
                monte_carlo += np.sum(experiment.data_without_mec)
                norm = np.sum(experiment.result) / monte_carlo
                diff = monte_carlo - np.sum(experiment.result)
                norm_term = self.norm_penalty * np.sum(experiment.result)
                norm_diff = diff / norm_term
            else:
                norm = 1
            print(experiment.name, experiment.covariance[0][0])
            unnormalized = calculate_score(matrix,
                                           experiment.monte_carlo * norm,
                                           experiment.result,
                                           experiment.error,
                                           experiment.covariance,
                                           experiment.data_without_mec * norm)
            if self.norm_penalty and 'minerva' in experiment.name:
                unnormalized[0] += norm_diff**2  # * experiment.bins
            if normalize_by_bins:
                score.append(np.array(unnormalized) / experiment.bins)
            else:
                score.append(np.array(unnormalized))

        return np.array(score)


@click.command()
@click.option('--iterations', default=20, type=int)
@click.option('--plot', is_flag=True)
@click.option('--show_matrix', is_flag=True)
@click.option('--population_size', default=10, type=int)
@click.option('--mutation_probability', default=0.05, type=float)
@click.option('--smoothness', default=0, type=float)
@click.option('--initial_matrix', default="ones", type=str)
@click.option('--eval_type', default="chi2_cov", type=str)
@click.option('--breed', is_flag=True)
@click.option('--breed_disable', default=0, type=int)
@click.option('--save', is_flag=True)
@click.option('--shape_covariance', is_flag=True)
@click.option('--normalize_by_cross_section', is_flag=True)
@click.option('--normalize_by_bins', is_flag=True)
@click.option('--pseudoinverse', is_flag=True)
@click.option('--normalization_penalty', default=0.0, type=float)
@click.option('--custom_covariance')
def combo_scaling(iterations: int,
                  plot: bool,
                  show_matrix: bool,
                  population_size: int,
                  mutation_probability: float,
                  initial_matrix: str,
                  smoothness: float,
                  eval_type: str,
                  breed: bool,
                  breed_disable: int,
                  save: bool,
                  shape_covariance: bool,
                  normalize_by_cross_section: bool,
                  normalize_by_bins: bool,
                  pseudoinverse: bool,
                  normalization_penalty: bool,
                  custom_covariance: str
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

    experiments = [('minerva_neutrino', 156), ('minerva_antyneutrino', 60),
                   ('T2K_neutrino', 58), ('T2K_antyneutrino', 58)]

    experimental_data = []
    for experiment, bins in experiments:
        data = utils.import_data_from_path(f"Data/{experiment}")
        experiment_info = utils.load_experimental_data(experiment)
        if "minerva" in experiment:
            events = utils.events_to_cross_section(data.events.reshape((24 * 24, -1)),
                                                   data.nof_events,
                                                   data.cross_section,
                                                   experiment_info[0].shape,
                                                   *experiment_info[4:6])
            if custom_covariance and experiment == 'minerva_neutrino':
                with open(custom_covariance, 'rb') as f:
                    covariance_matrix = np.linalg.pinv(pickle.load(f))

            if shape_covariance:
                if experiment == 'minerva_neutrino':
                    import covariance_matrix
                    covariance = covariance_matrix.covariance
                elif experiment == 'minerva_antyneutrino':
                    from minerva_antyneutrino import minerva_antyneutrino_covariance
                    covariance = minerva_antyneutrino_covariance

                if custom_covariance and experiment == 'minerva_neutrino':
                    with open(custom_covariance, 'rb') as f:
                        covariance = pickle.load(f)

                data = experiment_info[0].reshape(-1)
                shape_matrix, _, _ = utils.shape_matrix_extraction(data,
                                                                   covariance)
                covariance_matrix = shape_matrix

                if pseudoinverse:
                    covariance_matrix = np.linalg.pinv(covariance_matrix)
                else:
                    covariance_matrix = np.linalg.inv(covariance_matrix)
                experiment_info[2] = covariance_matrix
        else:
            events = data.events
        experimental_data.append(Experiment(experiment, events,
                                            experiment_info, bins))

    combined_score = Combined_score(experimental_data, normalization_penalty)

    # tmp_matrix = np.ones((24 * 24))

    score = combined_score(original_matrix, normalize_by_bins,
                           normalize_by_cross_section)

    print(f"Score without reweighting:\n{score}")

    population = [original_matrix]
    population += [mutate(original_matrix, mutation_probability, smoothness)
                   for _ in range(1, population_size)]

    if eval_type == 'chi2_cov':
        key = lambda m: np.sum(m[0][:, 0])
    elif eval_type == 'chi2_std':
        key = lambda m: np.sum(m[0][:, 1])
    elif eval_type == 'chi2_both':
        key = lambda m: np.sum(m[0])

    smoothness = smoothness if 0 < smoothness < 1 else 0

    start = time.time()
    try:
        for i_no in range(1, iterations + 1):
            chi2 = [combined_score(population[i], normalize_by_bins,
                                   normalize_by_cross_section)
                    for i in range(population_size)]

            population_with_score = list(zip(chi2, population))

            population_with_score.sort(key=key,
                                       reverse=False)

            print(i_no, population_with_score[0][0][:, 0],
                  sum(population_with_score[0][0][:, 0]))

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
            'iterations': i_no,
            'population_size': population_size,
            'mutation_probability': mutation_probability,
            'initial_matrix': initial_matrix,
            'smoothness': smoothness if smoothness else 1,
            'eval_type': eval_type,
            'breed': breed,
            'chi2_cov': save_chi2[:, 0],
            'chi2_std': save_chi2[:, 1]
        }
        summed_chi2 = np.sum(save_chi2, axis=0)
        name = f'scale_cov{summed_chi2[0]:.2f}_std{summed_chi2[1]:.2f}.pkl'
        with open(os.path.join('Combined_scaling', name), 'wb') as file:
            pickle.dump(population_with_score[0][1].reshape(24, 24), file)
            pickle.dump(config, file)
            print(f"Scaling matrix saved to {name}")


if __name__ == '__main__':
    combo_scaling()
