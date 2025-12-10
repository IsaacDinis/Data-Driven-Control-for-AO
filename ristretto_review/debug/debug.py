import numpy as np
import astropy.io.fits as fits
from matplotlib import pyplot as plt
from dd_utils import *
commands = fits.getdata("commands.fits")
res = fits.getdata("res.fits")
pol = fits.getdata("pol.fits")

commands_2 = fits.getdata("commands_vib.fits")
res_2 = fits.getdata("res_vib.fits")
pol_2 = fits.getdata("pol_vib.fits")
n_fft = 1000
fs = 2000

delay = 1

mode = 0
start = 200

pol_3 = pol_reconstruct(commands, -res, delay)
pol_4 = pol_reconstruct(commands_2, -res_2, delay)

# pol_3 = commands - res

# pol_4 = commands_2 - res_2

pol_psd, f, _ = compute_fft_mag_welch(pol[start:,mode], n_fft, fs)

pol_psd_2, f, _ = compute_fft_mag_welch(pol_2[start:,mode], n_fft, fs)


pol_psd_3, f, _ = compute_fft_mag_welch(pol_3[start:,mode], n_fft, fs)

pol_psd_4, f, _ = compute_fft_mag_welch(pol_4[start:,mode], n_fft, fs)

# plt.figure()
# plt.plot(commands[:,mode])
# plt.plot(res[:,mode])
# plt.plot(pol[:,mode])
# plt.legend(["com","res","pol"])



# plt.figure()
# plt.plot(commands[:,mode])
# plt.plot(commands_2[:,mode])

# plt.show()



plt.figure()
plt.plot(pol[:,mode])
plt.plot(pol_2[:,mode])

plt.figure()
plt.loglog(f,pol_psd)
plt.loglog(f,pol_psd_2)



plt.figure()
plt.loglog(f,pol_psd_3)
plt.loglog(f,pol_psd_4)

plt.show()
