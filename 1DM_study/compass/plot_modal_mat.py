import numpy as np
import matplotlib.pyplot as plt
M_DM0_2_M_DM1 = np.load('calib_mat/M_DM0_2_M_DM1.npy')

plt.figure()
plt.imshow(M_DM0_2_M_DM1)
plt.show()