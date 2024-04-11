import numpy as np
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
import astropy.io.fits as pfits

tilt_vibration = pfits.getdata('tilt_vibration.fits')

m_size = 404 #saxo+
s_size = 400

pad_size = int((m_size-s_size)/2)

pupil_grid = make_pupil_grid(s_size)
zernike_basis = make_zernike_basis(3, 1, pupil_grid)
piston = np.array(zernike_basis[0].shaped)
tilt = np.array(zernike_basis[1].shaped)
print(np.std(tilt,where = piston ==1))
plt.figure()
plt.imshow(tilt)

tilt_pad = np.pad(tilt,pad_size)
num_repeats = tilt_vibration.shape[0]
tilt_2D = np.dstack([tilt_pad])
pfits.writeto('tilt_2D.fits', tilt_2D, overwrite = True)
# tilt_sequence*= tilt_vibration[np.newaxis,np.newaxis,:]