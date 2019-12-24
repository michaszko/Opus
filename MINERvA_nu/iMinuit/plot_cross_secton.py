import pickle
from utilis import plot_matrix, plot_squares
import numpy as np
import matplotlib as plt
from minerva import matrix, cross_section

with open('../data/shape_matrix.pkl', 'rb') as f:
    shape_matrix = pickle.load(f)
    
with open('../data/minerva_neutrino_non_multiplied.pkl', 'rb') as f:
    minerva_neutrino_non_multiplied = pickle.load(f)
    
with open('../data/minerva_antyneutrino_non_multiplied.pkl', 'rb') as f:
    minerva_antyneutrino_non_multiplied = pickle.load(f)
    
with open('../data/T2K_neutrino_non_multiplied.pkl', 'rb') as f:
    T2K_neutrino_non_multiplied = pickle.load(f)
    
with open('../data/T2K_antyneutrino_non_multiplied.pkl', 'rb') as f:
    T2K_antyneutrino_non_multiplied = pickle.load(f)
    
with open('../data/minerva_neutrino_multiplied.pkl', 'rb') as f:
    minerva_neutrino_multiplied = pickle.load(f)
    
with open('../data/minerva_antyneutrino_multiplied.pkl', 'rb') as f:
    minerva_antyneutrino_multiplied = pickle.load(f)
    
with open('../data/T2K_neutrino_multiplied.pkl', 'rb') as f:
    T2K_neutrino_multiplied = pickle.load(f)
    
with open('../data/T2K_antyneutrino_multiplied.pkl', 'rb') as f:
    T2K_antyneutrino_multiplied = pickle.load(f)

experiments = [ minerva_neutrino_non_multiplied, minerva_neutrino_multiplied, 
                minerva_antyneutrino_non_multiplied, minerva_antyneutrino_multiplied,
                T2K_neutrino_non_multiplied, T2K_neutrino_multiplied,
                T2K_antyneutrino_non_multiplied,T2K_antyneutrino_multiplied]

for x in experiments:
  x[np.isnan(x)] = 0

min_nu_max = max(max(map(max, minerva_neutrino_non_multiplied)),max(map(max, minerva_neutrino_multiplied)))
min_anu_max = max(max(map(max, minerva_antyneutrino_non_multiplied)),max(map(max, minerva_antyneutrino_multiplied)))
t2k_nu_max = max(max(map(max, T2K_neutrino_non_multiplied)),max(map(max, T2K_neutrino_multiplied)))
t2k_anu_max = max(max(map(max, T2K_antyneutrino_non_multiplied)),max(map(max, T2K_antyneutrino_multiplied)))

plot_matrix(minerva_neutrino_non_multiplied,
            name='minerva_neutrino_non_multiplied',
            leg='MINER$\\nu$A $\\nu$ \nbefore',
            max_=min_nu_max)

plot_matrix(minerva_neutrino_multiplied,
            name='minerva_neutrino_multiplied',
            leg='MINER$\\nu$A $\\nu$ \nafter',
            max_=min_nu_max)

# Test -- cross-checking my and Tomek result on Minerva
# plot_matrix(np.matmul(np.ones((1,156)),matrix), 
#             name="minerva_nu", 
#             leg="MINER$\\nu$A $\\nu$")

# matrix[matrix!=matrix]=0
# print(np.sum(matrix)*2500/cross_section/1e6)

plot_matrix(minerva_antyneutrino_non_multiplied,
            name='minerva_antyneutrino_non_multiplied',
            leg='MINER$\\nu$A $\\bar{\\nu}$ \nbefore',
            max_=min_anu_max)

plot_matrix(minerva_antyneutrino_multiplied,
            name='minerva_antyneutrino_multiplied',
            leg='MINER$\\nu$A $\\bar{\\nu}$ \nafter',
            max_=min_anu_max)

plot_matrix(T2K_neutrino_non_multiplied,
            name='T2K_neutrino_non_multiplied',
            leg='T2K $\\nu$ \nbefore',
            max_=t2k_nu_max)

plot_matrix(T2K_neutrino_multiplied,
            name='T2K_neutrino_multiplied',
            leg='T2K $\\nu$ \nafter',
            max_=t2k_nu_max)

plot_matrix(T2K_antyneutrino_non_multiplied,
            name='T2K_antyneutrino_non_multiplied',
            leg='T2K $\\bar{\\nu}$ \nbefore',
            max_=t2k_anu_max)

plot_matrix(T2K_antyneutrino_multiplied,
            name='T2K_antyneutrino_multiplied',
            leg='T2K $\\bar{\\nu}$ \nafter',
            max_=t2k_anu_max)

plot_squares(shape_matrix, name="wq", leg="value 1")
