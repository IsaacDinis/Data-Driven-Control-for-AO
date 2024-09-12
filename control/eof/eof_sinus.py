import astropy.io.fits as pfits
import numpy as np
from matplotlib import pyplot as plt
import numpy.matlib
n_modes = 1
l = 10000
delay = 2
fs = 1000
max_freq = 50
freq = np.random.rand(n_modes,1)*max_freq*2*np.pi
T= l/fs
order = 40
t = np.arange(0,(order+l+delay)/fs,1/fs)
# data = np.sin(np.matlib.repmat(t,n_modes,1)*freq)
# data = pfits.getdata('P.fits').reshape(1,-1)*1000
data = 10*np.sin(np.matlib.repmat(t,n_modes,1)*2*np.pi*50+np.random.rand(1)*180)
plt.figure()
plt.plot(data.T)
plt.show()
D = np.zeros((order*n_modes,l))
P = np.zeros((n_modes,l))
for i in range(order):
    D[i*n_modes:(i+1)*n_modes,:] = data[:,delay+i:delay+i+l]
P = data[:,0:l]
D_inv = np.linalg.pinv(D.T)
F = (D_inv @ P.T).T
plt.figure()
plt.plot(F[0,:])
plt.show()
# data =
