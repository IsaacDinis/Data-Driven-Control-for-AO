import dao
import numpy as np

from dd4compass import K_dd

modes_shm =dao.shm('/tmp/papyrus_modes.im.shm')
dm = dao.shm('/tmp/dmCmd02.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()
n_modes = 195
M2V = M2V[:,:n_modes]
old_modes = np.zeros((M2V.shape[1],1),np.float32)
g = 0.6
K_dd = K_dd(5, 2, np.diag(np.ones(M2V.shape[1])), M2V, 1000000, 2000)

while True:
	modes = modes_shm.get_data(check = True)[:n_modes].squeeze()
	# command = g*modes + old_modes
	# voltage = -M2V@command
	voltage = K_dd.step(modes)
	dm.set_data(voltage.astype(np.float32))
	# old_modes = command