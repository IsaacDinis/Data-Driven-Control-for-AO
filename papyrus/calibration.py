import dao
import numpy as np
import time
from matplotlib import pyplot as plt
from dd4compass import K_dd


modes_shm =dao.shm('/tmp/papyrus_modes.im.shm')
slopes_shm = dao.shm('/tmp/papyrus_slopes.im.shm')
dm = dao.shm('/tmp/dmCmd02.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()

slopes = slopes_shm.get_data()
n_slopes = slopes.shape[0]

n_modes = 100
amp = 0.5

M2S = np.zeros((n_slopes,n_modes))

for i in range(n_modes):
    dm.set_data((M2V[:,i]*amp).astype(np.float32))
    time.sleep(0.1)
    slopes = slopes_shm.get_data()/amp
    M2S[:,i] += slopes.squeeze().copy()/2
    dm.set_data((-M2V[:,i]*amp).astype(np.float32))
    time.sleep(0.1)
    slopes = slopes_shm.get_data()/amp
    M2S[:,i] -= slopes.squeeze().copy()/2

S2M = np.linalg.pinv(M2S)
dao.shm('/tmp/S2M.shm',S2M.astype(np.float32))

amp = 0.5
mode_n = 1
dm.set_data((M2V[:,mode_n]*amp).astype(np.float32))
time.sleep(0.01)
slopes = slopes_shm.get_data(check = True, semNb = 5).squeeze()
modes = S2M@slopes
print(modes[mode_n]/amp)
plt.figure()
plt.plot(modes/amp)
plt.show()