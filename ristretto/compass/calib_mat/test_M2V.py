import numpy as np
from matplotlib import pyplot as plt
import astropy.io.fits as pfits

M2V = pfits.getdata('M2V_DM1.fits') # mode to volts matrix [n_act x n_modes]
plt.figure()
plt.plot(np.std(M2V,axis = 1))  
plt.show()

plt.figure()
plt.plot(M2V[0,:])  
plt.plot(M2V[10,:])  
plt.show()