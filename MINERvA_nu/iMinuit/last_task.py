import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.special     import expit
from scipy.stats       import chi2
from minerva import *

with open('../data/shape_matrix_old.pkl', 'rb') as f:
    shape_matrix = pickle.load(f)
    
covariance_noninverted_short = covariance_noninverted*1e80
shape_matrix_short = shape_matrix*1e80
cross_daniel_short = np.reshape(np.transpose(np.reshape(cross_daniel,(12,13))),(1,156))*1e40
cross_daniel_good = np.reshape(np.transpose(np.reshape(cross_daniel,(12,13))),(1,156))[0]

for i in np.flip(np.where(cross_daniel_good==0)[0],axis=0):
    shape_matrix_short = np.delete(shape_matrix_short, (i), axis=0)
    shape_matrix_short = np.delete(shape_matrix_short, (i), axis=1)
    covariance_noninverted_short = np.delete(covariance_noninverted_short, (i), axis=0)
    covariance_noninverted_short = np.delete(covariance_noninverted_short, (i), axis=1)
    cross_daniel_short = np.delete(cross_daniel_short, i)
    
shape_matrix_short_inv = np.linalg.pinv(shape_matrix_short)
shape_matrix_inv = np.linalg.pinv(shape_matrix)

# savetxt("../data/daniel_covariance_short.dat", covariance_noninverted_short)
# savetxt("../data/cross_daniel_short.dat", cross_daniel_short)
# savetxt("../data/shape_matrix_short_inv.dat", shape_matrix_short_inv)

sample = 10000
df = 156
zeros = len(np.where(cross_daniel_good==0)[0])
eff_df = df-zeros-1

random_vector = np.random.multivariate_normal(cross_daniel_good*1e40, covariance_noninverted*1e80, sample)
chi2_vector = np.random.chisquare(eff_df,sample)
other_vector = np.random.poisson(eff_df,size=sample)
random_vector_norm = random_vector/np.resize(np.sum(random_vector,axis=1),(sample,1))*np.sum(cross_daniel_good*1e40)

diff = (cross_daniel_good*1e40 - random_vector_norm)
result = np.diag(np.matmul(np.matmul(diff,shape_matrix_inv*1e-80),np.transpose(diff)))

fig, ax = plt.subplots(1, 1)

weights = np.ones_like(result)/len(result)

range1, range2, bins, density = (70, 250, 30, 1)

# ax.hist(other_vector, histtype="stepfilled", density=density, weights=weights,
#         label="poisson with " + str(df-zeros-1) + " dof", range=(range1, range2), bins=bins, alpha=0.6)

ax.hist(chi2_vector, histtype="stepfilled", density=density, weights=weights,
        label="$\\chi2$ with " + str(df-zeros-1) + " dof", range=(range1, range2), bins=bins, alpha=0.6)

ax.hist(result, histtype="stepfilled", density=density, weights=weights,
        label="$\\chi2$ of results", range=(range1, range2), bins=bins, alpha=0.6)

ax.axvline(df-zeros-1, color='r', label=str(df-zeros-1))

ax.set_xlabel("$\\chi$2")
ax.set_ylabel("normalized")
ax.set_title("Shape_matrix")
plt.legend()

plt.show()
# fig.savefig("shape_matrix_norm_chi2.pdf", dpi=70, bbox_inches='tight', pad_inches=0, quality=95)