import time
import dao
import numpy as np




modes_shm =dao.shm('/tmp/papyrus_modes.im.shm')
dm = dao.shm('/tmp/dmCmd02.im.shm')
dmTurb=dao.shm('/tmp/dmCmd03.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()
slopes_shm = dao.shm('/tmp/papyrus_slopes.im.shm')
S2M = dao.shm("/tmp/S2M.shm").get_data()
n_modes = 1
M2V = M2V[:,:n_modes]


command = np.zeros((n_modes,1),np.float32)
command[0] = 0.1
voltage = M2V@command

eps = 1e-3
dm.set_data(voltage.astype(np.float32))

while (S2M@slopes_shm.get_data(check = True, semNb = 5))[0] < command[0] - eps:
    dmTurb.get_data(check = True, semNb = 5)
command[0] = -0.1
voltage = M2V@command
dm.set_data(voltage.astype(np.float32))

time_start = time.time()
while (S2M@slopes_shm.get_data(check = True, semNb = 5))[0] > command[0] + eps:
    dmTurb.get_data(check = True, semNb = 5)
latency = time.time() - time_start
command[0] = 0
voltage = M2V@command
dm.set_data(voltage.astype(np.float32))


print(latency)


