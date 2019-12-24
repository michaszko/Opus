import os
import numba
import pickle
import numpy as np
import importlib
from typing import Tuple, List
from types import ModuleType
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D


def import_data_from_path(data_path: str) -> ModuleType:
    """Returns module `events` located at given path"""
    module_path = str(os.path.join(data_path, "events")).replace("/", ".")

    return importlib.import_module(module_path)


def load_experimental_data(experiment: str
                           ) -> List[np.ndarray]:
    with open(f"{experiment}.pkl", 'rb') as f:
        result = pickle.load(f)
        err = pickle.load(f)
        covariance = pickle.load(f)
        without_mec = pickle.load(f)
        trans_left = pickle.load(f)
        long_left = pickle.load(f)

    return [result, err, covariance, without_mec, trans_left, long_left]


@numba.jit
def events_to_cross_section(events: np.ndarray,
                            nof_events: int,
                            cross_section: float,
                            target_shape: Tuple[int, int],
                            trans_left: List[int],
                            long_left: List[int],
                            bin_size_term: bool=True
                            ) -> np.ndarray:
    """Calculates cross sections in each bin of given histogramized events"""
    events_ = events.reshape((24, 24, target_shape[0], target_shape[1]))
    cross_sections = np.zeros((24, 24, target_shape[0], target_shape[1]))
    norm = cross_section / nof_events * 1e6 * 12/13
    for long in range(target_shape[1]):
        for trans in range(target_shape[0]):
            bin_size = (long_left[long + 1] - long_left[long])
            bin_size *= (trans_left[trans + 1] - trans_left[trans])
            bin_size = bin_size if bin_size_term else 1
            for q in range(24):
                for w in range(24):
                    # sigma = events_[q][w][trans][long] * norm / bin_size
                    sigma = events_[q][w][trans][long] * norm / 2500
                    cross_sections[q][w][trans][long] = sigma
    return cross_sections.reshape((24 * 24, target_shape[0] * target_shape[1]))


@numba.jit
def chi2(experiment: np.ndarray,
         nuwro: np.ndarray,
         covariance_matrix: np.ndarray,
         ) -> float:
    """Calculate chi2 error between experimental data and Monte Carlo"""
    chi2 = 0

    for j in range(experiment.shape[0]):
        for k in range(experiment.shape[0]):
            err = experiment[j] - nuwro[j]
            err *= covariance_matrix[j][k]
            err *= experiment[k] - nuwro[k]
            chi2 += err

    return chi2


@numba.jit
def chi2_standard(experiment: np.ndarray,
                  nuwro: np.ndarray,
                  target_err: np.ndarray
                  ) -> float:
    chi2 = 0

    for j in range(experiment.shape[0]):
        if target_err[j]:
            err = (experiment[j] - nuwro[j])
            err /= target_err[j]
            chi2 += err ** 2

    return chi2


def shape_matrix_extraction(data: np.ndarray,
                            covariance: np.ndarray):
    nT = np.sum(data)
    cov_row = np.sum(covariance, axis=1)
    cov_sum = np.sum(cov_row)

    shape_mat = np.zeros_like(covariance)
    norm_mat = np.ones_like(covariance) * cov_sum / nT**2

    for i in range(len(data)):
        for j in range(len(data)):
            norm_mat[i][j] *= data[i] * data[j]
            shape_mat[i][j] = covariance[i][j] + norm_mat[i][j]
            shape_mat[i][j] -= (data[i] * cov_row[i] +
                                data[j] * cov_row[j]) / nT

    mixed_mat = covariance - shape_mat - norm_mat

    return shape_mat, mixed_mat, norm_mat


def show_matrix(matrix: np.ndarray, logscale: bool=False) -> None:
    from matplotlib import pyplot

    pyplot.rcParams['font.family'] = 'serif'
    # pyplot.rcParams['font.sans-serif'] = 'Tahoma'
    pyplot.rcParams['font.size'] = 18

    m = matrix.copy()
    for i in range(24):
        for j in range(i + 1, 24):
            m[j][i] = None

    fig, ax = pyplot.subplots()
    if logscale:
        from matplotlib.colors import LogNorm
        vmin = np.min(matrix)
        vmax = np.max(matrix)
        ax.matshow(m, origin='lower',
                   norm=LogNorm(vmin=vmin, vmax=vmax))
    else:
        ax.matshow(m, origin='lower',
                   vmax=2, vmin=0.0)

    for (i, j), z in np.ndenumerate(m):
        if not np.isnan(z):
            ax.text(j, i, '{:0.1f}'.format(z), ha='center', va='center')

    pyplot.show()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--show_matrix', type=str,
                        help="Path to pickled matrix")
    parser.add_argument('--logscale', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--data_location')
    parser.add_argument('--experiment', default="minerva_neutrino")
    parser.add_argument('--save', action='store_true')
    parser.add_argument('--write_to_C')
    args = parser.parse_args()

    if args.write_to_C:
        import pickle
        with open(args.write_to_C, 'rb') as f:
            matrix = pickle.load(f)
        with open('scaling_matrix.h', 'w') as f:
            f.write("double scaling_matrix[24][24] = {\n\t")
            for row in matrix:
                for i in row:
                    f.write(f'{i}, ')
            f.write('\n}\n')

    if args.show_matrix:
        from matplotlib import pyplot
        pyplot.rc('text', usetex=True)
        import pickle
        with open(args.show_matrix, 'rb') as matrix_file:
            matrix = pickle.load(matrix_file)
            try:
                matrix_info = pickle.load(matrix_file)
                print(matrix_info)
            except:
                pass

        show_matrix(matrix, args.logscale)

        if args.plot:
            if not args.data_location:
                data_path = os.path.join("Data", args.experiment)
                data = import_data_from_path(data_path)
            else:
                data = import_data_from_path(args.data_location)

            experiment_info = load_experimental_data(args.experiment)
            target_result, target_err = experiment_info[:2]
            covariance_matrix, data_without_mec = experiment_info[2:4]
            trans_left, long_left = experiment_info[4:6]

            target_shape = target_result.shape
            print(target_shape)

            events = data.events.reshape((24 * 24, -1))
            if "minerva" in args.experiment:
                print('MINERvA data')
                events = events_to_cross_section(events,
                                                 data.nof_events,
                                                 data.cross_section,
                                                 target_shape,
                                                 trans_left,
                                                 long_left)
            else:
                data_without_mec *= 12 / 13
                data.events = data.events * 12 / 13

            pyplot.matshow(np.sum(events, axis=-1).reshape(24, 24),
                           origin='lower')
            pyplot.colorbar()
            pyplot.show()

            target_shape = target_result.shape
            scaled_data = np.dot(matrix.reshape(-1),
                                 events).reshape(target_shape)

            pyplot.matshow(np.sum(events, axis=-1).reshape(24, 24) * matrix,
                           origin='lower')
            # pyplot.title(f'Scaled MINERvA nu alpha=0.3')
            pyplot.colorbar()
            pyplot.show()

            before = np.sum(events)
            after = np.sum(np.sum(events, axis=-1).reshape(24, 24) * matrix)
            non_mec = np.sum(data_without_mec)
            print(f"before: {before}, after: {after}, change: {after/before - 1}")
            print(f"data: {np.sum(target_result)}, non MEC: {non_mec}")
            print(f"diff before: {1-(before+non_mec)/np.sum(target_result)}, after: {1-(after+non_mec)/np.sum(target_result)}")
            # chi2_score = chi2(target_result,
            #                   scaled_data + data_without_mec,
            #                   covariance_matrix)
            # print(f"chi2: {chi2_score}")
            # non_scaled = events_to_cross_section(np.dot(np.ones(576),
            #                                             events
            #                                             ).reshape(target_shape),
            #                                      data.nof_events,
            #                                      data.cross_section,
            #                                      target_shape,
            #                                      trans_left,
            #                                      long_left)

            non_scaled = np.dot(np.ones(576), events).reshape(target_shape)
            particle = "$\\bar{\\nu}_{\\mu}$" if 'anty' in args.experiment else "$\\nu_{\\mu}$"

            if "T2K" in args.experiment:
                sklejka = [
                    [14],
                    [1, 1, 1, 1, 10],
                    [1, 1, 1, 1, 2, 8],
                    [1, 1, 1, 1, 2, 8],
                    [1, 1, 1, 1, 2, 2, 6],
                    [1, 1, 1, 1, 2, 2, 2, 4],
                    [2, 1, 1, 2, 3, 2, 3],
                    [2, 1, 1, 2, 2, 1, 1, 1, 1, 2],
                    [3, 2, 2, 2, 2, 1, 1, 1]
                ]

                data_bin = 0
                for i, bins in enumerate(sklejka):
                    fig, ax = pyplot.subplots(1)
                    ax.clear()
                    start_bin = 0
                    target = []
                    error = []
                    n_scaled = []
                    scaled = []
                    rectangles = []
                    mid = []
                    width = []
                    edges = [5]
                    final_bin = data_bin + len(bins)
                    target.append(target_result[data_bin:final_bin])
                    error.append(target_err[data_bin:final_bin])
                    n_scaled.append(
                        non_scaled[data_bin:final_bin] + data_without_mec[data_bin:final_bin])
                    scaled.append(
                        scaled_data[data_bin:final_bin] + data_without_mec[data_bin:final_bin])
                    for m_bin in bins:
                        m = trans_left[start_bin] + \
                            trans_left[start_bin + m_bin]
                        mid.append(m / 2)
                        w = -trans_left[start_bin] + \
                            trans_left[start_bin + m_bin]
                        width.append(w / 2)
                        start_bin += m_bin
                        edges.append(trans_left[start_bin])

                    for k,j in enumerate(range(data_bin, final_bin)):
                        rect = Rectangle((edges[k], 0),
                        2*width[k], data_without_mec[j])
                        rectangles.append(rect)

                    pc = PatchCollection(rectangles, color="black", alpha=0.1, edgecolor="none")
                    ax.add_collection(pc)

                    data_bin = final_bin

                    pyplot.errorbar(mid, [*target[0]],
                                    xerr=width, yerr=[*error[0]],
                                    fmt='k.', ecolor='k')
                    l1=pyplot.hlines([*target[0]],
                                  xmin=edges[:-1], xmax=edges[1:], colors='k')
                    l2=pyplot.hlines([*n_scaled[0]],
                                  xmin=edges[:-1], xmax=edges[1:], colors='b',
                                  linestyles='dashed')
                    l3=pyplot.hlines([*scaled[0]],
                                  xmin=edges[:-1], xmax=edges[1:], colors='r',
                                  linestyle='dotted')
                    print([*n_scaled[0]][:-1])
                    pyplot.vlines(edges[1:-1],
                                  ymin=[*n_scaled[0]][:-1], ymax=[*n_scaled[0]][1:], colors='b',
                                  linestyles='dashed')
                    pyplot.vlines(edges[1:-1],
                                  ymin=[*scaled[0]][:-1], ymax=[*scaled[0]][1:], colors='r',
                                  linestyle='dotted')
                    ax.autoscale()
                    l4=Line2D([0, 1], [0, 1], color="black", alpha=0.1, linewidth=10)
                    ax.legend([l1,l2,l3,l4],[f"T2K {particle}", "Default NuWro", 'Scaled NuWro','non-MEC'],
                                  loc='upper right', fontsize='medium')
                    pyplot.xlabel("$p$ [MeV]",fontsize='medium')
                    pyplot.ylabel("$\\sigma$ [cm$^{2}$/GeV$^{2}$/nucleon]", fontsize='medium')
                    pyplot.title(f"{long_left[i]} $<$ cos $<$ {long_left[i+1]}", fontsize='medium')
                    pyplot.xlim(left=0, right=3200) #powinno być 5200
                    if args.save:
                        pyplot.savefig(f"{i+1}.pdf", bbox_inches='tight', pad_inches=0)
                    pyplot.show()

                pyplot.show()
                import sys
                sys.exit(0)

            # target_result = target_result.T
            # target_err = target_err.T
            # data_without_mec = data_without_mec.T
            # target_shape = target_result.shape

            # trans_left = long_left

            mid = np.array(trans_left[:-1]) + np.array(trans_left[1:])
            mid = mid / 2
            width = np.array(trans_left[1:]) - np.array(trans_left[:-1])
            width = width / 2

            target_result = target_result.T
            data_without_mec = data_without_mec.T
            scaled_data = scaled_data.T
            non_scaled = non_scaled.T
            target_shape = target_result.shape
            if args.experiment == 'minerva_neutrino':
                target_err = target_err.T.reshape(target_shape)

            # print(target_result.shape)
            # print(scaled_data.shape)
            # print(data_without_mec.shape)

            # target_result = target_result.T

            if args.experiment == 'minerva_antyneutrino':
                target_err = target_err.reshape((6, 10))
                target_err = target_err.T

            for i in range(target_shape[0]):
                fig, ax = pyplot.subplots(1)
                rectangles = []

                for k,j in enumerate(range(len(data_without_mec[i]))):
                    rect = Rectangle((trans_left[k], 0),
                    2*width[k], data_without_mec[i][k], angle=0)
                    rectangles.append(rect)

                pc = PatchCollection(rectangles, color="black", alpha=0.1, edgecolor="none")

                ax.add_collection(pc)

                sc = scaled_data[i] + data_without_mec[i]
                nsc = non_scaled[i] + data_without_mec[i]
                # experimental uncertainty
                pyplot.errorbar(x=mid, y=target_result[i],
                                yerr=target_err[i], xerr=width,
                                fmt='k.', ecolor='k')
                # horizontal
                l1=pyplot.hlines(target_result[i],
                              xmin=trans_left[:-1],
                              xmax=trans_left[1:])
                l2=pyplot.hlines(nsc,
                              xmin=trans_left[:-1],
                              xmax=trans_left[1:], colors='b',
                              linestyles='dashed')
                l3=pyplot.hlines(sc,
                              xmin=trans_left[:-1],
                              xmax=trans_left[1:], colors='r',
                              linestyles='dotted')
                # vertical
                pyplot.vlines(trans_left[1:-1],
                              ymin=nsc[:-1], ymax=nsc[1:], colors='b',
                              linestyles='dashed')
                pyplot.vlines(trans_left[1:-1],
                              ymin=sc[:-1], ymax=sc[1:], colors='r',
                              linestyles='dotted')
                ax.autoscale()
                l4=Line2D([0, 1], [0, 1], color="black", alpha=0.1, linewidth=10)
                ax.legend([l1,l2,l3,l4],[f"MINERvA {particle}", "Default NuWro", 'Scaled NuWro','non-MEC'],
                              loc='upper right', fontsize='medium')
                pyplot.title(f"{long_left[i]} $< p_L$ [MeV] $<$ {long_left[i+1]}", fontsize='medium')
                pyplot.xlabel("$p_T$ [MeV]", fontsize='medium')
                pyplot.ylabel("$\\sigma$ [cm$^{2}$/GeV$^{2}$/nucleon]", fontsize='medium')
                if args.save:
                    pyplot.savefig(f"{i+1}.pdf", bbox_inches='tight', pad_inches=0)
                pyplot.show()
