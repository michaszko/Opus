import numpy as np

mec_dim       = 24  # dimension of matrix of rescaling
events_number = 5e6  # number of events in sample
cross_section = 1.19683e-39  # calculated cross section in .txt file
alpha         = 0.3 # parameter of smoothness

delta_pT           = np.array([75, 75, 100, 75, 75, 75, 75, 150, 150, 150, 250, 250, 1000])
delta_pT_plot      = np.array([75, 75, 100, 75, 75, 75, 75, 150, 150, 150, 250, 250, 1000, 0])
delta_pL           = np.array([500, 500, 500, 500, 500, 500, 500, 1000, 2000, 2000, 5000, 5000])
pT                 = np.array([75, 150, 250, 325, 400, 475, 550, 700, 850, 1000, 1250, 1500, 2500])
pT_plot            = np.array([75, 150, 250, 325, 400, 475, 550, 700, 850, 1000, 1250, 1500, 2500, 3000])
pL                 = np.array([1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 8000, 10000, 15000, 20000])
bins_width_vanilla = np.kron(delta_pT, delta_pL)
bins_width         = np.outer(bins_width_vanilla,  np.ones(mec_dim**2))

# loading files
matrix_vanilla      = np.fromfile("../data/matrix_LFG.dat",                sep=" ")  # distribution of MEC events in bins 
cross_nuwro         = np.fromfile("../data/cros_total_nuwro_LFG.dat",      sep=" ")  # cross section of nuwro (w-o MEC)
cross_daniel        = np.fromfile("../data/cros_total_daniel.dat",         sep=" ")  # cross section of MINERvA
covariance_vanilla  = np.fromfile("../data/daniel_covariance_vanilla.dat", sep=" ")  # covariance matrix
cross_error_vanilla = np.fromfile("../data/cros_error_daniel.dat",         sep=" ")  # error bars

# modifying files
matrix_unnormalized    = np.reshape(matrix_vanilla,(156, mec_dim**2))  # reshaping MEC matrix
matrix                 = matrix_unnormalized * cross_section * 1e6 * 12 / 13 / bins_width / events_number  # normalizing matrix

covariance_noninverted = np.reshape(covariance_vanilla, (156,156))  # reshaping covarince matrix
covariance             = np.linalg.pinv(covariance_noninverted)  # pseudo-inverse covariance matrix (used later); inverse doesnt exist

cross_error            = np.where(cross_error_vanilla==0, 1, cross_error_vanilla)  # removing zeros from error bars

# rotation of results due to covariance matrix's look 
cross_daniel_cov       = np.reshape(cross_daniel,(12,13)).transpose().flatten()  
cross_nuwro_cov        = np.reshape(cross_nuwro,(12,13)).transpose().flatten()