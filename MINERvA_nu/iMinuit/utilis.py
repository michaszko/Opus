'''Definitions of useful functions'''
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import RegularPolyCollection
from math import pi
from scipy.special import expit
from minerva import alpha, matrix, cross_daniel_cov, cross_nuwro_cov, covariance

def plot_matrix(x, name="lol", leg='', maxx=None, ylab=True, show=True):
    """Drawing matrix and saving it to .pdf file"""
    from matplotlib import rc
    rc('text', usetex=True)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 18

    rescale = np.reshape(x, (24, 24))

    for i in range(24):
        for j in range(i + 1, 24):
            rescale[j][i] = None

    fig, ax = plt.subplots()

    if maxx==None:
        im = ax.imshow(rescale, origin="lower", extent=(1, 1200, 1, 1200))   
    else:
        im = ax.imshow(rescale, origin="lower", extent=(1, 1200, 1, 1200), vmax=maxx)
    cbar = ax.figure.colorbar(im, ax=ax)
    # cbar.ax.set_title("$\\left[\\frac{\mathrm{cm}^2}{\mathrm{GeV}^2}\\right]$",
    #                     ha="left", x=1, fontsize='small')
    if ylab==True:
        cbar.ax.set_ylabel("$\\frac{\\mathrm{d}^2\\sigma}{\\mathrm{d}q\\mathrm{d}\\omega}$" +
                           "$\\left[\\frac{\\mathrm{cm}^2}{\\mathrm{GeV}^2}\\right]$")
    ax.set_xlabel("q [MeV]")
    ax.set_ylabel("$\\omega$ [MeV]")
    ax.set_xticks(np.arange(0, 1400, 200))
    ax.set_xticklabels(np.arange(0, 1201, 200))
    ax.annotate(leg, xy=(10, 255), xycoords='axes points',
                size=14, ha='left', va='top',
                bbox=dict(boxstyle='round', fc='w'))

    fig.savefig("plots/" + name + ".pdf", dpi=300,
                bbox_inches='tight', pad_inches=0.1, quality=95)
    if show==True:
        plt.show()

def plot_squares(data, name="wq", leg=''):
    """Function for plotting matrices with squares"""
    from matplotlib import rc
    rc('text', usetex=True)
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 18
    
    factor = 20
    pow_factor = 0.6
    M = np.array((data))**(pow_factor)*factor
    M[M > 80] = 80

    fig, ax = plt.subplots(1, 1, subplot_kw={'aspect': 'equal'})
    ax.set_xlim(-0.5, M.shape[0] - 0.5)
    ax.set_ylim(-0.5, M.shape[1] - 0.5)

    # xy locations of each square center
    xy = np.indices(M.shape)[::-1].reshape(2, -1).T
    w = (M.ravel())
    ec = RegularPolyCollection(4, sizes=w, rotation=pi/4, facecolor='black',
                               edgecolor='black', offsets=xy, 
                               transOffset=ax.transData, norm=10)

    ax.add_collection(ec)
    ax.set_xticks(np.arange(0, 28, 4))
    ax.set_xticklabels(np.arange(0, 1201, 200))
    ax.set_yticks(np.arange(0, 28, 4))
    ax.set_yticklabels(np.arange(0, 1201, 200))
    ax.set_xlabel("q [MeV]")
    ax.set_ylabel("$\\omega$ [MeV]")
    plt.legend((leg,), loc=2, markerscale=.9)

    fig.savefig("plots/" + name + ".pdf", dpi=300,
                bbox_inches='tight', pad_inches=0.1, quality=95)
    plt.show()