import numpy as np
from matplotlib import pyplot as plt

V2V = np.load('calib_mat/V_DM0_2_V_DM1.npy')
n_act_LODM = V2V.shape[1]
command_LODM = np.zeros(n_act_LODM)
command_LODM[int(n_act_LODM/2)] = 1
command_HODM = V2V@command_LODM
command_HODM[command_HODM<0.01] = 0
plt.plot(command_HODM)
plt.show()