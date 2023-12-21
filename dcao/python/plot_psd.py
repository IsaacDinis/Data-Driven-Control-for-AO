import numpy as np
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
import compute_psd
import astropy.io.fits as pfits
import utils
from scipy import signal

modal_res = pfits.getdata('../dcao_ol/saxoplus_KL_res.fits')
mode = 0
fs = 2760
n_average = 60
window_size = 200
# psd, freq = compute_psd.compute_psd_fft(modal_res, n_average, window_size,fs)
psd, freq, psd_map = compute_psd.compute_psd_welch(modal_res, window_size,fs)
# psd /= np.max(psd)
psd_log = 10*np.log10(psd)
fig, ax = plt.subplots(constrained_layout=True)

ax.set_title('mode {:d}'.format(mode))
ax.set_xscale('log')
ax.set_ylabel("mag. [dB]")
ax.set_xlabel("freq. [Hz]")
ax.plot(freq,psd_log[:,mode])


# fig1, ax1= plt.subplots(constrained_layout=True)
#
# ax1.set_title('mode {:d}'.format(mode))
# ax1.set_ylabel("freq. [Hz]")
# ax1.set_xlabel("time")
# # ax1.set_yscale('log')
#
# ax1.imshow(psd_map[:,:,0])
# #
plt.show()

pfits.writeto('psd.fits', psd, overwrite = True)
pfits.writeto('freq.fits', freq, overwrite = True)