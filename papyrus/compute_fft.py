import dao
import numpy as np
import time
from matplotlib import pyplot as plt
from dd4compass import K_dd
from dd_utils import *


modes_shm =dao.shm('/tmp/papyrus_modes.im.shm')
dm = dao.shm('/tmp/dmCmd02.im.shm')
dm_turb_shm=dao.shm('/tmp/dmCmd03.im.shm')
res_buf_shm = dao.shm('/tmp/res_buf.shm')
command_buf_shm = dao.shm('/tmp/command_buf.shm')
pol_buf_shm = dao.shm('/tmp/pol_buf.shm')
turb_buf_shm = dao.shm('/tmp/turb_buf.shm')


M2V = dao.shm("/tmp/m2c.im.shm").get_data()
buf_size = 1024
n_modes = 1
n_fft = 300

turb_buf_fft = np.zeros((int(n_fft/2),n_modes),np.float32)
pol_buf_fft = np.zeros((int(n_fft/2),n_modes),np.float32)
command_buf_fft = np.zeros((int(n_fft/2),n_modes),np.float32)
res_buf_fft = np.zeros((int(n_fft/2),n_modes),np.float32)
f_buff = np.zeros((int(n_fft/2),1),np.float32)

turb_buf_fft_shm = dao.shm('/tmp/turb_buf_fft.shm',turb_buf_fft.astype(np.float32))
pol_buf_fft_shm = dao.shm('/tmp/pol_buf_fft.shm',pol_buf_fft.astype(np.float32))
command_buf_fft_shm = dao.shm('/tmp/command_buf_fft.shm',command_buf_fft.astype(np.float32))
res_buf_fft_shm = dao.shm('/tmp/res_buf_fft.shm',res_buf_fft.astype(np.float32))
f_buff_shm =  dao.shm('/tmp/f_buff.shm',f_buff.astype(np.float32))

M2V = M2V[:,:n_modes]
V2M = np.linalg.pinv(M2V)


fs = 100

latency = 1
old_time = time.time()



while True:
    time.sleep(5)
    pol_buf = pol_buf_shm.get_data()
    turb_buf = turb_buf_shm.get_data()
    res_buf = res_buf_shm.get_data()
    command_buf = command_buf_shm.get_data()

    pol_fft, f, _ = compute_fft_mag_welch(pol_buf, n_fft, fs)
    turb_fft, f, _ = compute_fft_mag_welch(turb_buf, n_fft, fs)
    res_fft, f, _ = compute_fft_mag_welch(res_buf, n_fft, fs)
    command_fft, f, _ = compute_fft_mag_welch(command_buf, n_fft, fs)

    if n_modes == 1:
        pol_fft = pol_fft[:,np.newaxis]
        turb_fft = turb_fft[:,np.newaxis]
        res_fft = res_fft[:,np.newaxis]
        command_fft = command_fft[:,np.newaxis]
        f = f[:,np.newaxis]
    # turb_fft = np.ones((int(n_fft/2)+1,n_modes))
    # pol_buf_fft_shm.set_data(np.log10(pol_fft[:-1]).astype(np.float32))
    # turb_buf_fft_shm.set_data(np.log10(turb_fft[:-1]).astype(np.float32))
    # res_buf_fft_shm.set_data(np.log10(res_fft[:-1]).astype(np.float32))
    # command_buf_fft_shm.set_data(np.log10(command_fft[:-1]).astype(np.float32))
    # f_buff_shm.set_data(f[:-1].astype(np.float32))

    pol_buf_fft_shm.set_data((pol_fft[1:,:]).astype(np.float32))
    turb_buf_fft_shm.set_data((turb_fft[1:,:]).astype(np.float32))
    res_buf_fft_shm.set_data((res_fft[1:,:]).astype(np.float32))
    command_buf_fft_shm.set_data((command_fft[1:,:]).astype(np.float32))
    f_buff_shm.set_data(f[1:,:].astype(np.float32))