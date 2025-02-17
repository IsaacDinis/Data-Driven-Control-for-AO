import numpy as np
import astropy.io.fits as pfits
from matplotlib import pyplot as plt
from hcipy.field import make_pupil_grid
from hcipy.mode_basis import make_zernike_basis
from scipy.interpolate import interpn

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
plt.scatter(act_pos[648,0],act_pos[648,1],marker='.', color="blue")
plt.ylim((0,pupil_diam))
plt.xlim((0,pupil_diam))
plt.ylim(max(plt.ylim()), min(plt.ylim()))
plt.show()

#######################################################################################################
n_act_cut = 43
pupil_cut_diam = act_size*n_act_cut
act_center_cut = 1092

act_pos_cut = np.squeeze(act_pos[np.hypot((act_pos[:,0]-act_pos[act_center_cut,0]),(act_pos[:,1]-act_pos[act_center_cut,1]))<330].copy())
act_index_cut = np.squeeze(act_index[np.hypot((act_pos[:,0]-act_pos[act_center_cut,0]),(act_pos[:,1]-act_pos[act_center_cut,1]))<330].copy())

plt.figure()
plt.imshow(flat)
plt.scatter(act_pos_cut[:,0],act_pos_cut[:,1],marker='.', color="red")
plt.scatter(act_pos_cut[375,0],act_pos_cut[375,1],marker='.', color="blue")
plt.ylim((0,pupil_diam))
plt.xlim((0,pupil_diam))
plt.ylim(max(plt.ylim()), min(plt.ylim()))
plt.show()

##############################################################################################################

dummy = np.ceil(act_pos[act_center_cut,0]-pupil_cut_diam/2)

flat_cut = flat[int(np.ceil(act_pos[act_center_cut,1]-pupil_cut_diam/2)):int(np.ceil(act_pos[act_center_cut,1]+pupil_cut_diam/2)),
           int(np.ceil(act_pos[act_center_cut,0]-pupil_cut_diam/2)):int(np.ceil(act_pos[act_center_cut,0]+pupil_cut_diam/2))].copy()

pupil_grid_cut = make_pupil_grid(np.ceil(pupil_cut_diam), 1)
zernike_basis_cut = make_zernike_basis(1, 1, pupil_grid_cut)
pupil_valid_cut = zernike_basis_cut[0].shaped
act_pos_cut2 = np.array((act_pos_cut[:,0]-np.ceil(act_pos[act_center_cut,0]-pupil_cut_diam/2),act_pos_cut[:,1]-np.ceil(act_pos[act_center_cut,1]-pupil_cut_diam/2))).T
# flat_cut *= pupil_valid_cut
plt.figure()
plt.imshow(flat_cut)
plt.scatter(act_pos_cut2[:,0],act_pos_cut2[:,1],marker='.', color="red")
plt.show()
##############################################################################################
m_size = 336
s_size = 320
pad_size = int((m_size-s_size)/2)


flat_size = flat_cut.shape[0]
dummy_x,dummy_y = np.meshgrid(np.linspace(0,flat_size-1,s_size),np.linspace(0,flat_size-1,s_size))


flat_interp = interpn((np.arange(flat_size), np.arange(flat_size)), flat_cut,(dummy_x, dummy_y ))

flat_pad = np.pad(flat_interp,pad_size)
flat_record = np.dstack([flat_pad])
flat_record *= 1e6

plt.figure()
plt.imshow(flat_pad)
plt.show()
pfits.writeto('flat.fits',flat_record, overwrite=True)
pfits.writeto('act_pos.fits',act_pos_cut, overwrite=True)

m_size = 124

flat_size = flat.shape[0]
flat_array = np.arange(flat_size)
dummy_x,dummy_y = np.meshgrid(np.linspace(0,flat_size-1,m_size),np.linspace(0,flat_size-1,m_size))
flat_interp = interpn((flat_array, flat_array), flat,(dummy_x, dummy_y ))
flat_interp *= 1e6
flat_piston = np.mean(flat_interp[flat_interp!=0])
flat_interp -= flat_piston
flat_record = np.dstack([flat_interp])
pfits.writeto('flat.fits',flat_record, overwrite=True)

plt.figure()
plt.imshow(flat)

plt.figure()
plt.imshow(flat_interp)
