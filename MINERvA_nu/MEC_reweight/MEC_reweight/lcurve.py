import os
import pickle
import numpy as np
from matplotlib import pyplot
from matplotlib.colors import LogNorm

from utils import show_matrix

alpha = []
chi2 = []
chi2_std = []

# pyplot.rc('text', usetex=True)
# pyplot.xkcd()

for matrix in os.listdir():
    try:
        with open(matrix, 'rb') as file:
            m = pickle.load(file)
            vmin = np.min(m)
            vmax = np.max(m)
            info = pickle.load(file)

            relative = 0
            allowed = 0
            for i in range(24):
                for j in range(i + 1):
                    relative += m[j][i]
                    allowed += 1

            for i in range(24):
                for j in range(i + 1, 24):
                    m[j][i] = None
            # pyplot.matshow(m, origin='lower',
            #                norm=LogNorm(vmin=vmin, vmax=vmax))
            # print(vmin, vmax, relative / allowed)
            pyplot.imshow(m, origin='lower', extent=(0, 1200, 0, 1200))
            # pyplot.matshow(m, origin='lower')
            pyplot.colorbar()
            pyplot.title(f"$\\alpha: {info['smoothness']}$")
            pyplot.xticks(np.linspace(0, 1200, 7))
            pyplot.yticks(np.linspace(0, 1200, 7))
            pyplot.xlabel('$q$')
            pyplot.ylabel('$\\omega$')
            pyplot.show()
        # if info['experiment'] == "minerva_neutrino":
        # print(info)
        chi2.append(np.sum(info['chi2_cov']))
        chi2_std.append(np.sum(info['chi2_std']))
        alpha.append(info['smoothness'] if info['smoothness'] else 1)
    except:
        pass


pyplot.show()
pyplot.show()
pyplot.show()

# # alpha: [0.1, 0.3, 0.5, 0.2, 1, 0.4]
# non_norm_chi2_std = [317.40502181939524, 284.64886644063364, 490.5365344399719]
# non_norm_chi2_cov = [261.2818080081728, 247.35950536576257,
#                      294.6034527922023, 337.3261636645773, 382.6093896872222]
# pyplot.plot(non_norm_chi2_std, [0.1, 0.2, 0.3], 'y.')
# pyplot.plot(non_norm_chi2_cov, [0.1, 0.2, 0.3, 0.4, 0.5], 'g.')

# alpha.append(0)
# chi2.append(309.86679306329415)
# chi2_std.append(383.97918650927767)


# titles = {
#     "minerva_antyneutrino": "MINERvA $\\bar{\\nu}_{\\mu}$",
#     "minerva_neutrino": "MINERvA $\\nu_{\\mu}$",
#     "T2K_antyneutrino": "T2K $\\bar{\\nu}_{\\mu}$",
#     "T2K_neutrino": "T2K $\\nu_{\\mu}$",
#     "Combined_scaling": "All experiments",
# }

pyplot.plot(alpha, chi2, '.', markersize=15)
# pyplot.plot(alpha, chi2_std, 'r.')
# pyplot.legend(['with covariance', "w/o covariance"])
pyplot.title("$\\chi^2(\\alpha)$ dependency")
pyplot.xlabel("$\\alpha$")
pyplot.ylabel("$\\chi^2$")
# pyplot.imsave('lcurve.pdf')
pyplot.savefig('lcurve.pdf')
pyplot.show()
