import dao
import numpy as np
import time
from matplotlib import pyplot as plt
from dd4compass import K_dd

modes_shm =dao.shm('/tmp/papyrus_modes.im.shm')
dm = dao.shm('/tmp/dmCmd02.im.shm')
dmTurb=dao.shm('/tmp/dmCmd03.im.shm')
M2V = dao.shm("/tmp/m2c.im.shm").get_data()

buf_size = 1024
n_modes = 1
res_buf = np.zeros((buf_size,n_modes),np.float32)
command_buf = np.zeros((buf_size,n_modes),np.float32)
res_buf_shm = dao.shm('/tmp/res_buf.shm',res_buf.astype(np.float32))
command_buf_shm = dao.shm('/tmp/command_buf.shm',command_buf.astype(np.float32))


M2V = M2V[:,:n_modes]
old_command = np.zeros((M2V.shape[1],1),np.float32)
g = 0.2
training_size = 3000
fs = 100
# K_dd = K_dd(5, 1, np.diag(np.ones(M2V.shape[1])), M2V, training_size, fs)
new_time = time.time()
while True:
    # old_time = new_time
    # new_time = time.time()
    # print(new_time-old_time)
    dmTurb.get_data(check = True, semNb = 5)
    res_buf = res_buf_shm.get_data()
    command_buf = command_buf_shm.get_data()
    res_buf = np.roll(res_buf, -1, axis=0)
    command_buf = np.roll(command_buf, -1, axis=0)
    # time.sleep(0.001)
    # modes = modes_shm.get_data(check = True)[:n_modes]
    modes = modes_shm.get_data(check = True, semNb = 5)[:n_modes].squeeze()
    command = g*modes + old_command
    # command*= 0
    voltage = -M2V@command
    command_buf[-1, :] = command
    res_buf[-1, :] = modes
    res_buf_shm.set_data(res_buf.astype(np.float32))
    command_buf_shm.set_data(command_buf.astype(np.float32))
    
    # voltage = K_dd.step(modes)
    dm.set_data(voltage.astype(np.float32))

    old_command = command