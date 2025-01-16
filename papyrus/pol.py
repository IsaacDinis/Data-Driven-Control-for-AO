import dao
import numpy as np
import time
from matplotlib import pyplot as plt
from dd4compass import K_dd
from dynamic_array_plot import SingleArrayDynamicPlot, MultiArrayDynamicPlot

def pol_reconstruct(command_buff, measurement_buff, delay):
    delay_floor = int(np.floor(delay))
    delay_ceil = int(np.ceil(delay))
    delay_frac,_ = np.modf(delay)
    if delay_ceil == delay_floor:
        pol = measurement_buff[-1,:] + command_buff[-delay_ceil-1,:]
        # print(command_buff[-delay_ceil-1,0])
        # print(measurement_buff[-1, 0])
    else:
        pol = measurement_buff[-1, :] + (1 - delay_frac) * command_buff[-delay_floor-1, :] + delay_frac * command_buff[-delay_ceil-1,:]
    return pol


modes_shm =dao.shm('/tmp/papyrus_modes.im.shm')
dm = dao.shm('/tmp/dmCmd02.im.shm')
dm_turb_shm=dao.shm('/tmp/dmCmd03.im.shm')
res_buf_shm = dao.shm('/tmp/res_buf.shm')
command_buf_shm = dao.shm('/tmp/command_buf.shm')

M2V = dao.shm("/tmp/m2c.im.shm").get_data()
buf_size = 1024
n_modes = 1
turb_buf = np.zeros((buf_size,n_modes),np.float32)
pol_buf = np.zeros((buf_size,n_modes),np.float32)
turb_buf_shm = dao.shm('/tmp/turb_buf.shm',turb_buf.astype(np.float32))
pol_buf_shm = dao.shm('/tmp/pol_buf.shm',pol_buf.astype(np.float32))

M2V = M2V[:,:n_modes]
V2M = np.linalg.pinv(M2V)


fs = 100

latency = 1
old_time = time.time()

while True:
    dm_turb = dm_turb_shm.get_data(check = True, semNb = 6)
    pol_buf = pol_buf_shm.get_data()
    turb_buf = turb_buf_shm.get_data()
    res_buf = res_buf_shm.get_data()
    command_buf = command_buf_shm.get_data()
    pol_buf = np.roll(pol_buf, -1, axis=0)
    turb_buf = np.roll(turb_buf, -1, axis=0)
    turb = V2M@dm_turb
    pol = pol_reconstruct(command_buf, res_buf, latency)
    pol_buf[-1, :] = pol
    turb_buf[-1, :] = turb
    pol_buf_shm.set_data(pol_buf.astype(np.float32))
    turb_buf_shm.set_data(turb_buf.astype(np.float32))
    if(time.time()-old_time > 30):
        old_time = time.time()
        print(np.std(res_buf)/np.std(pol_buf))