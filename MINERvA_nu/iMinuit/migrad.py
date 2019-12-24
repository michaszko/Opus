import numpy as np
from iminuit import Minuit
from minerva import *
import random
from scipy.special import expit
from utilis import plot_matrix
import pickle


def mysigmo(x, trans=0):
    '''my modified sigmoid'''
    return 100 * expit(x * 200 - trans)


def starting_matrix(x, rand=True):
    '''starting matrix - 1 above diagonal, x below diagonal'''
    tmp = []
    if x == -1:
        tmp = pickle.load(open("matrices/matrix_366.0.pkl", "rb"))
    else:
        if rand == True:  # with random fluctuation
            for i in range(0, 24):
                tmp += [1] * i + [
                    x + 0.1 * random.random() for a in range(24 - i)
                ]
        else:  # without random fluctuation
            for i in range(0, 24):
                tmp += [1] * i + [x] * (24 - i)

    return tmp


def punishment(rescale):
    upper = []  # elements' indexes over diagonal
    for i in range(24):
        for j in range(0, i):
            upper += [24 * i + j]

    lower = []  # elements' indexes below diagonal
    for i in range(24):
        for j in range(i, 24):
            lower += [24 * i + j]

    # d - diagonal_line r - right_line b - bottom line i - inner . - both
    #
    #     ---------
    #     |      .|
    #     |    dir|
    #     |  diiir|
    #     |.bbbbb.|
    #     ---------

    bottom_line = range(24)
    right_line = range(23, 24**2, 24)
    diagonal_line = [24 * i + i for i in range(24)]
    inner = [
        x for x in lower if x not in bottom_line and x not in right_line
        and x not in diagonal_line
    ]

    pun = 0  # punishment

    for i in bottom_line[1:-1]:  # for bottom side of matrix
        pun += mysigmo(-rescale[i] /
                       max(rescale[i + 1], rescale[i + 24], rescale[i - 1]) +
                       (1 - alpha))
        pun += mysigmo(rescale[i] /
                       min(rescale[i + 1], rescale[i + 24], rescale[i - 1]) -
                       1 / (1 - alpha))
    for i in right_line[1:-1]:  # for right side of matrix
        pun += mysigmo(-rescale[i] /
                       max(rescale[i + 24], rescale[i - 1], rescale[i - 24]) +
                       (1 - alpha))
        pun += mysigmo(rescale[i] /
                       min(rescale[i + 24], rescale[i - 1], rescale[i - 24]) -
                       1 / (1 - alpha))
    for i in diagonal_line[1:-1]:  # for diagonal side of matrix
        pun += mysigmo(-rescale[i] / max(rescale[i + 1], rescale[i - 24]) +
                       (1 - alpha))
        pun += mysigmo(rescale[i] / min(rescale[i + 1], rescale[i - 24]) - 1 /
                       (1 - alpha))
    for i in inner[0:-1]:  # for the rest
        pun += mysigmo(-rescale[i] / max(rescale[i + 1], rescale[
            i - 1], rescale[i + 24], rescale[i - 24]) + (1 - alpha))
        pun += mysigmo(rescale[i] / min(rescale[i + 1], rescale[
            i - 1], rescale[i + 24], rescale[i - 24]) - 1 / (1 - alpha))

    return pun


def chi2_cov(rescaling, pun=True):
    '''definition of chi^2 with covariance'''
    mec = np.reshape(np.matmul(matrix, rescaling),
                     (12, 13)).flatten()  # rescaled MEC

    # chi^2 = tmp.cov.tmp
    tmp = (cross_daniel_cov - cross_nuwro_cov - mec)
    if pun == True:
        return (np.matmul(np.matmul(tmp, covariance), tmp) +
                punishment(rescaling))
    else:
        return (np.matmul(np.matmul(tmp, covariance), tmp))


if __name__ == '__main__':
    # every element should not be negative
    non_negative = np.reshape((0.1, None) * mec_dim**2, (mec_dim**2, 2))

    # fixing over diagonal elements
    triangular = []
    for i in range(0, 24):
        triangular += [True] * i + [False] * (24 - i)

    m = Minuit.from_array_func(
        chi2_cov,
        starting_matrix(-1),
        error=0.1,
        errordef=1,
        limit=non_negative,
        # fix=triangular,
        use_array_call=True)

    m.migrad(ncall=10000)
    plot_matrix(m.np_values(),
                "migrad_" + 
                str(round(chi2_cov(m.np_values(), pun=False), 2)),
                str(round(chi2_cov(m.np_values(), pun=False), 2)),
                ylab=False,
                show=False)
    pickle.dump(
        m.np_values(),
        open(
            "matrices/matrix_" + str(round(chi2_cov(m.np_values()), 2)) +
            ".pkl", "wb"))

    for i in range(100):
        m = Minuit.from_array_func(
            chi2_cov,
            m.np_values(),
            error=0.1,
            errordef=1,
            limit=non_negative,
            # fix=triangular,
            use_array_call=True)

        m.migrad(ncall=10000)
        plot_matrix(m.np_values(),
                    "migrad_" +
                    str(round(chi2_cov(m.np_values(), pun=False), 2)),
                    str(round(chi2_cov(m.np_values(), pun=False), 2)),
                    ylab=False,
                    show=False)
        pickle.dump(
            m.np_values(),
            open(
                "matrices/matrix_" + str(round(chi2_cov(m.np_values()), 2)) +
                ".pkl", "wb"))
