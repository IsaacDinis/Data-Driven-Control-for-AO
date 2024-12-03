import subprocess
import signal
import os
import time
from scipy.io import loadmat
import dao
import numpy as np
import astropy.io.fits as pfits

time.sleep(1) # wait for SHM to be created
amp=10
turb = pfits.getdata("tilt_traj.fits")
dmTurb=dao.shm('/tmp/dmCmd03.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()
# Infinite loop to wait for Ctrl+C
fs = 100
while True:  
    for k in range(turb.shape[0]):
        dmTurb.set_data(M2V[:,0]*turb[k].astype(np.float32)*amp)
        time.sleep(1/fs)
    for k in np.linspace(turb.shape[0]-1,0,turb.shape[0]):
        dmTurb.set_data(M2V[:,0]*turb[int(k)].astype(np.float32)*amp)
        time.sleep(1/fs)
