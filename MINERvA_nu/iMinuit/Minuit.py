
# get_ipython().run_line_magic('pylab', 'inline')
from pprint import pprint # we use this to pretty print some stuff later

import numpy as np
import pandas as pd

from matplotlib.colors      import LogNorm
from matplotlib             import pyplot as plt
from matplotlib.collections import RegularPolyCollection

from iminuit import Minuit
from iminuit import minimize

# from scipy.optimize import minimize
from scipy.special     import expit
from scipy.stats       import chi2
from scipy.interpolate import make_interp_spline, BSpline

from  minerva import *

# # checking Minuit -- easier method

# res = minimize(chi2_cov, 
#                starting_matrix(1), 
#                bounds=bnd,
#                options={"maxfev":10000000})

# plot_matrix(res["x"], chi2_cov(res["x"], False) , "migrad", res["success"])

# for i in range(4):
#     res = minimize(chi2_cov, 
#                    res["x"], 
#                    bounds=bnd,
#                    options={"maxfev":10000000})

#     plot_matrix(res["x"], chi2_cov(res["x"], False) , "migrad", res["success"])


# In[ ]:


# # how does the result change with respect to starting matrix
# chi2_of_starting_param = []
# starting_param = []

# for i in numpy.arange(0, 1, 0.1):
#     tmp = m(i)
#     tmp.migrad()
#     chi2_of_starting_param += [tmp.fval]
#     starting_param += [i]
    
# xlabel("param of starting matrix")
# ylabel("$\\chi^2$")
# plot(starting_param, chi2_of_starting_param,".")
# # savefig("/home/michaszko/Desktop/lol.pdf")


# In[ ]:


# # how does the result osccilate 
# chi2_of_starting_param = []
# starting_param = []

# for i in range(10):
#     tmp = m(0.1)
#     tmp.migrad()
#     chi2_of_starting_param += [tmp.fval]
#     starting_param += [i]
    
# xlabel("param of starting matrix")
# ylabel("$\\chi^2$")
# plot(starting_param, chi2_of_starting_param,".")
# # savefig("/home/michaszko/Desktop/lol.pdf")


# # minimize() 
# ### with scipy
# bounds and constrains
tmp = []
for i in range(0,24):  # fixing over diagonal elements
    tmp += [1,1] * i + [0.01, None] * (24-i)
    
bnd = np.reshape(tmp,(576,2))  # bounds - no zero in order to avoid division by it

##################################

upper = []
for i in range(24):
    for j in range(0, i):
        upper += [24*i + j]  # elements' indexes over diagonal 

lower = []
for i in range(24):
    for j in range(i, 24):
        lower += [24*i + j]  # elements' indexes below diagonal 
        
# d - diagonal_line 
# r - right_line
# b - bottom line
# i - inner
# . - both 
#
#     ---------
#     |      .|
#     |    dir|
#     |  diiir|
#     |.bbbbb.|
#     ---------

bottom_line    = range(24)
right_line     = range(23, 24 ** 2, 24)
diagonal_line  = [24*i + i for i in range(24)]
inner          = [x for x in lower if x not in bottom_line and x not in right_line and x not in diagonal_line]


alpha = 0.3  # parameter of smoothness 

cons = []

# for i in upper:  # for upper elements of matrix
#     cons += [{'type': 'eq', 'fun': lambda x: x[i] - 1}]

for i in bottom_line[1:-1]:  # for bottom side of matrix
    cons += [{'type': 'ineq', 'fun': lambda x: x[i] / max( \
            x[i + 1], \
            x[i + 24],\
            x[i - 1]) - (1 - alpha)}]
    cons += [{'type': 'ineq', 'fun': lambda x:  -x[i] / min( \
            x[i + 1], \
            x[i + 24],\
            x[i - 1]) + 1/(1 - alpha) }]
    
for i in right_line[1:-1]:  # for right side of matrix
    cons += [{'type': 'ineq', 'fun': lambda x: x[i] / max( \
            x[i + 24], \
            x[i - 1],\
            x[i - 24]) - (1 - alpha)}]
    cons += [{'type': 'ineq', 'fun': lambda x:  -x[i] / min( \
            x[i + 24], \
            x[i - 1],\
            x[i - 24]) + 1/(1 - alpha) }]
    
for i in diagonal_line[1:-1]:  # for diagonal side of matrix
    cons += [{'type': 'ineq', 'fun': lambda x: x[i] / max( \
            x[i + 1], \
            x[i - 24]) - (1 - alpha)}]
    cons += [{'type': 'ineq', 'fun': lambda x:  -x[i] / min( \
            x[i + 1], \
            x[i - 24]) + 1/(1 - alpha) }]
    
for i in inner[0:-1]:  # for the rest 
    cons += [{'type': 'ineq', 'fun': lambda x: x[i] / max( \
            x[i + 1],\
            x[i - 1],\
            x[i + 24],\
            x[i - 24]) - (1 - alpha)}]
    cons += [{'type': 'ineq', 'fun': lambda x:  -x[i] / min( \
            x[i + 1],\
            x[i - 1],\
            x[i + 24],\
            x[i - 24]) + 1/(1 - alpha) }]# possible methods
# L-BFGS-B, TNC
# COBYLA, SLSQP, trust-constr

mth = "trust-constr"

# rescale_matrix = starting_matrix(0.05)

res = minimize(chi2_cov, 
#                    rescale_matrix, 
               starting_matrix(1),
               method=mth, 
               bounds=bnd, 
               constraints=cons, 
               options={"maxiter":2000})

#     rescale_matrix = res["x"]

savetxt("longer_1.dat", res["x"])

plot_matrix(res["x"], res["fun"])# ploting results from file
ress = np.array(loadtxt("Minuit/matrices/migrad_punish_0.3.dat")).flatten()
plot_matrix(ress, chi2_cov(ress),titleTrue=False)
# plot_results(ress, "rescale_5")# # plotting results from recent calculations
# plot_matrix(res["x"], res["fun"])
# plot_results(res["x"], "rescale")# # how the results differ with changes of starting matrix?
# chi2_of_starting_param = []
# starting_param = []

# for j in numpy.arange(0, 1, 0.2):
#     res = minimize(chi2_cov, starting_matrix(j), method='SLSQP', bounds=bnd, constraints=cons)
#     chi2_of_starting_param += [res['fun']]
#     starting_param += [j]# # plot of upper results 
# plot(starting_param, chi2_of_starting_param,".")
# xlabel("param of starting matrix")
# ylabel("$\\chi^2$")
# savefig("/home/michaszko/Desktop/lol.pdf")
# xlim(0,2)



# # Plotting charts



# # Last task

# In[98]:




