import subprocess
import signal
import os
import time
from scipy.io import loadmat
import dao
import numpy as np
import astropy.io.fits as pfits


amp=50
turb = pfits.getdata("tilt_traj.fits")
dmTurb=dao.shm('/tmp/dmCmd03.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()
slopes_shm = dao.shm('/tmp/papyrus_slopes.im.shm')
# Infinite loop to wait for Ctrl+C
fs = 100
while True:  
    for k in range(turb.shape[0]):
        # full_turb = turb[k]+np.sin(2*pi*4)
        dmTurb.set_data(M2V[:,0]*turb[k].astype(np.float32)*amp)
        slopes_shm.get_data(check=True, semNb=10)
    for k in np.linspace(turb.shape[0]-1,0,turb.shape[0]):
        dmTurb.set_data(M2V[:,0]*turb[int(k)].astype(np.float32)*amp)
        slopes_shm.get_data(check=True, semNb=10)
