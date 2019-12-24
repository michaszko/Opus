import os
import pickle
from utils import show_matrix

folders = [
    # "Minerva_antyneutrino",
    # "Minerva_neutrino",
    # "T2K_neutrino",
    # "T2K_antyneutrino"
    "optimize_test"
]


# for f in folders:
for file in os.listdir(f"Matrices/normalized_minerva"):
    if '.pkl' in file:
        with open(f"Matrices/normalized_minerva/{file}", 'rb') as m:
            matrix = pickle.load(m)
            try:
                info = pickle.load(m)
                if info['experiment'] == 'minerva_neutrino':
                    print(f"Matrices/normalized_minerva/{file}", flush=True)
                    # show_matrix(matrix)
                    # for x in matrix.reshape(-1):
                    #     print(f'{x}, ', end='')
            except:
                continue
