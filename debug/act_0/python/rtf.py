import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
import compute_psd

ol = pfits.getdata('../compass/results/ol/zernike_res.pfits')
cl = pfits.getdata('../compass/results/cl/zernike_res.pfits')
n_average = 10
window_size = 390
fs = 2000

psd_average_ol, f = compute_psd.compute_psd_xcorr(ol,n_average,window_size,fs)
psd_average_ol = psd_average_ol[1:,:]
psd_average_cl, f = compute_psd.compute_psd_xcorr(cl,n_average,window_size,fs)
psd_average_cl = psd_average_cl[1:,:]

f = f[1:]
plt.semilogx(f,10*np.log10(psd_average_ol[:,0]))
plt.semilogx(f,10*np.log10(psd_average_cl[:,0]))
plt.semilogx(f,10*np.log10(psd_average_cl[:,0])-10*np.log10(psd_average_ol[:,0]))
plt.show()
