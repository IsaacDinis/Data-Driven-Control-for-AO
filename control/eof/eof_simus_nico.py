import astropy.io.fits as pfits
import numpy as np
from matplotlib import pyplot as plt
import numpy.matlib
from scipy.linalg import svd, pinv
n_modes = 1
l = 10000
delay = 2
fs = 1000
max_freq = 50
freq = np.random.rand(n_modes,1)*max_freq*2*np.pi
order = 50

data = pfits.getdata('P.fits').reshape(1,-1)

plt.figure()
plt.plot(data.T)
plt.show()

D = np.zeros((order*n_modes,l))
P = np.zeros((n_modes,l))


for i in range(l):
    D[:,i] = data[0,i:i+order]
P = data[0,order+delay:order+delay+l]



# Perform SVD (with econ='full_matrices=False' for economy-size decomposition)
u, s, v = svd(D, full_matrices=False)  # u * s * v.T will give back D



# Normalize singular values
snorm = s / np.max(s)
noise =  1
# Set threshold based on noise
threshold = 1e-4 if noise == 0 else 1e-3

threshold = 0
# Find nmax (index of minimum value close to threshold)
nmax = np.argmin(np.abs(snorm - threshold))

# Create sc matrix similar to s, using zero matrix and filling it with s
sc = np.zeros_like(s) + s
# sc[nmax:] = np.inf  # Commented based on the MATLAB code

# Compute Dp
Dp = v[:nmax, :].T @ pinv(np.diag(sc[:nmax])) @ u[:, :nmax].T

# Uncomment if you want Tikohonov regularization instead:
# lambda_ = 0.1  # Define lambda if needed
# Dp = np.linalg.inv(D.T @ D + lambda_ * np.eye(k)) @ D.T

# Compute Fi
F = P @ Dp

# D_inv = np.linalg.pinv(D.T)
# F = (D_inv @ P.T).T


plt.figure()
plt.plot(F)
plt.show()

plt.figure()
plt.plot(P.squeeze())
plt.show()
# data =