import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid
from hcipy.mode_basis import make_zernike_basis

pupil_diam = 752
pupil_grid = make_pupil_grid(pupil_diam, 1)
zernike_basis = make_zernike_basis(1, 1, pupil_grid)
pupil_valid = zernike_basis[0].shaped


flat = np.loadtxt('flat.txt')
pix_size = 26.63e-6 # pix/m
act_size_m = 400e-6 # m
act_size = act_size_m/pix_size
dm_size = np.floor(20e-3/pix_size)
n_act = 50

flat_size_0 = flat.shape[0]
flat_size_1 = flat.shape[1]
dummy_0 = int((flat_size_0 - dm_size) / 2)
dummy_1 = int((flat_size_1 - dm_size) / 2)
flat = flat[dummy_0:-dummy_0, dummy_1:-dummy_1]
flat *= pupil_valid
dummy= np.arange(0,dm_size,act_size)

dummy = np.linspace(act_size/2,pupil_diam-act_size/2,n_act)

dummy = np.meshgrid(dummy,dummy)
act_pos = np.array([dummy[0].flatten(),dummy[1].flatten()]).T


act_pos = act_pos[np.hypot((act_pos[:,0]-pupil_diam/2),(act_pos[:,1]-pupil_diam/2))<385] #385 to get 2040 actuators
act_index = np.arange(act_pos.shape[0])
plt.figure()
plt.imshow(flat)
plt.scatter(act_pos[:,0],act_pos[:,1],marker='.', color="red")
# plt.scatter(act_pos[1092,0],act_pos[1092,1],marker='.', color="blue")
plt.ylim((0,pupil_diam))
plt.xlim((0,pupil_diam))
plt.ylim(max(plt.ylim()), min(plt.ylim()))
plt.show()

pfits.writeto("act_pos.fits", act_pos, overwrite = True)

