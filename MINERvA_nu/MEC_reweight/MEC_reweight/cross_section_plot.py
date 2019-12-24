import os
import numba
import pickle
import importlib
import numpy as np
from matplotlib import pyplot

from utils import import_data_from_path, load_experimental_data, events_to_cross_section

# pyplot.xkcd()
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--show_matrix', type=str,
                        help="Path to pickled matrix")
    args = parser.parse_args()

    with open(args.show_matrix, 'rb') as matrix_file:
        matrix = pickle.load(matrix_file)
        try:
            matrix_info = pickle.load(matrix_file)
            print(matrix_info)
        except:
            pass

    plot, axes = pyplot.subplots(nrows=4)
    experiments = ['minerva_neutrino', 'minerva_antyneutrino',
                   'T2K_neutrino', 'T2K_antyneutrino']
    titles = ['MINERvA $\\nu_{\\mu}$', 'MINERvA $\\bar{\\nu}_{\\mu}$',
              'T2K $\\nu_{\\mu}$', 'T2K $\\bar{\\nu}_{\\mu}$']
    for i, (experiment, title) in enumerate(zip(experiments, titles)):
        data_path = os.path.join("Data", experiment)
        data = import_data_from_path(data_path)

        experiment_info = load_experimental_data(experiment)
        target_result, target_err = experiment_info[:2]
        covariance_matrix, data_without_mec = experiment_info[2:4]
        trans_left, long_left = experiment_info[4:6]

        target_shape = target_result.shape
        print(target_shape)

        events = data.events.reshape((24 * 24, -1))
        if "minerva" in experiment:
            print('MINERvA data')
            events = events_to_cross_section(events,
                                             data.nof_events,
                                             data.cross_section,
                                             target_shape,
                                             trans_left,
                                             long_left)

        x = axes[i].imshow(np.multiply(matrix,np.sum(events, axis=-1).reshape(24, 24)), # without rescaling
        # x = axes[i].imshow(np.sum(events, axis=-1).reshape(24, 24),					# with rescaling
                           origin='lower', extent=(0, 1200, 0, 1200))
                           # vmax=8.8399518421e-40)
        axes[i].set_title(title)
        axes[i].set_xlabel('$q$')
        axes[i].set_ylabel('$\\omega$')
        axes[i].set_xticks(np.linspace(0, 1200, 5))
        axes[i].set_yticks(np.linspace(0, 1200, 5))
        plot.colorbar(x, ax=axes[i])
        # axes[i].colorbar()

        output = open(experiments[i]+'_non_multiplied'+'.pkl', 'wb')
        pickle.dump(np.sum(events, axis=-1).reshape(24, 24), output)                        #without rescaling
        output.close()      
        output = open(experiments[i]+'_multiplied'+'.pkl', 'wb')
        pickle.dump(np.multiply(matrix,np.sum(events, axis=-1).reshape(24, 24)), output)     #with rescaling 
        output.close()

    pyplot.show()
    # pyplot.savefig('test.pdf')
