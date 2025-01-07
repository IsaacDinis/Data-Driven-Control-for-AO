import numpy as np
import astropy.io.fits as pfits
from dd_utils import *

amp=10

turb = pfits.getdata("tilt_traj.fits")

plt.figure()
plt.plot(turb)


n_fft = 1000
fs = 100
turb_fft,f,_ = compute_fft_mag_welch(turb, n_fft, fs)
plt.figure()
plt.loglog(f,turb_fft)

plt.show()
