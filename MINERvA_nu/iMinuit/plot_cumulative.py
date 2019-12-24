import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
from matplotlib import pyplot as plt

xaxis = np.array([250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000])
numu0 = np.array([4.72083e-41, 3.49173e-40, 7.67323e-40, 1.0775e-39 , 1.16089e-39, 1.17817e-39, 1.18393e-39, 
                  1.19822e-39, 1.20163e-39, 1.20835e-39, 1.2157e-39 , 1.21844e-39, 1.22084e-39, 1.22062e-39,
                  1.22703e-39, 1.22631e-39, 1.23066e-39, 1.23418e-39, 1.23726e-39, 1.23525e-39])
numu1 = np.array([2.5904e-41, 2.17321e-40 , 1.27653e-39  , 1.91771e-39  ,  2.03751e-39 ,  2.05693e-39 ,  
                  2.08052e-39 ,  2.10781e-39 ,  2.11098e-39 ,  2.11314e-39 ,  2.12311e-39 ,   2.12964e-39,
                  2.13434e-39,   2.14523e-39,   2.14614e-39,   2.15001e-39,   2.15064e-39,   2.15549e-39, 
                  2.16015e-39,   2.16047e-39])
anumu0 = np.array([1.80243e-41  ,1.09107e-40  ,2.3892e-40 ,3.87475e-40  , 5.13212e-40 , 6.13392e-40 , 
                   6.8724e-40  , 7.49625e-40 , 7.99818e-40 , 8.44939e-40 , 8.82372e-40 ,  9.10561e-40,  
                   9.38343e-40,  9.61103e-40,  9.81321e-40,  9.9705e-40 ,  1.01276e-39,  1.02822e-39,  
                   1.04152e-39,  1.05251e-39])
anumu1 = np.array([1.03895e-41, 6.46271e-41, 2.50655e-40, 5.27047e-40,  7.6247e-40,  9.53049e-40, 
                   1.10068e-39,  1.21156e-39,  1.31182e-39,  1.39162e-39,  1.4598e-39,   1.51675e-39,
                   1.57223e-39,   1.61816e-39,   1.65554e-39,  1.67746e-39,   1.71123e-39,   1.74383e-39,
                   1.76804e-39,   1.79152e-39])

xnew = np.linspace(xaxis.min(), xaxis.max(), 300)
xnew_1 = np.linspace(300, xaxis.max(), 300)
spl_1 = make_interp_spline(xaxis, numu0, k=1)
power_smooth_1 = spl_1(xnew)
spl_2 = make_interp_spline(xaxis, numu1, k=1)
power_smooth_2 = spl_2(xnew)
spl_3 = make_interp_spline(xaxis, anumu0, k=1)
power_smooth_3 = spl_3(xnew)
spl_4 = make_interp_spline(xaxis, anumu1, k=1)
power_smooth_4 = spl_4(xnew)


plt.plot(xnew,power_smooth_1, color="black", linestyle="dashed")
plt.plot(xnew,power_smooth_2, color="black", linestyle="default")
plt.xlabel("Neutrino energy [MeV]")
plt.ylabel("Cross section (per nucleon)")
plt.xlim(0)
plt.ylim(0)
plt.legend(('Valencia model', 'Phenomenological model'),loc=4, fontsize='medium')
plt.show()
plt.savefig("plots/" + "neutrino" + ".pdf", dpi=300, bbox_inches='tight', pad_inches=0.1, quality=95)


plt.plot(xnew,power_smooth_3, color="black", linestyle="dashed")
plt.plot(xnew,power_smooth_4, color="black", linestyle="default")
plt.xlabel("Antineutrino energy [MeV]")
plt.ylabel("Cross section (per nucleon)")
plt.xlim(0)
plt.ylim(0)
plt.legend(('Valencia model', 'Phenomenological model'),loc=4, fontsize='medium')
plt.show()
plt.savefig("plots/" + "antineutrino" + ".pdf", dpi=300, bbox_inches='tight', pad_inches=0.1, quality=95)
