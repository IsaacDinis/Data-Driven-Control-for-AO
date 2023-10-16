import numpy as np
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
import compute_psd
import astropy.io.fits as pfits
import utils
from scipy import signal

modal_res = pfits.getdata('../dcao_ol/zernike_saxoplus_res.fits')
mode = 0
fs = 2000
n_average = 10
window_size = 200
psd, freq = compute_psd.compute_psd_xcorr(modal_res, n_average, window_size,fs)
# psd /= np.max(psd)
psd = 10*np.log10(psd)
fig, ax = plt.subplots(constrained_layout=True)

ax.set_title('mode {:d}'.format(mode))
ax.set_xscale('log')
ax.set_ylabel("mag. [dB]")
ax.set_xlabel("freq. [Hz]")
ax.plot(freq,psd[:,mode])



f, Pxx_den = signal.periodogram(modal_res[:,0], fs)

plt.plot(f[1:], 10*np.log10(Pxx_den[1:]))

# plt.ylim([1e-3, 1e4])

plt.xlabel('frequency [Hz]')

plt.ylabel('PSD [V**2/Hz]')

plt.show()