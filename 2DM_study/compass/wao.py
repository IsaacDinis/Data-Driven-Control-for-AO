M2V_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM1.npy')
M2V_DM0 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM1.npy')
M2V = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V.npy')
M_DM0_2_M_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M_DM0_2_M_DM1.npy')
nact_DM0 = 88
V2V = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM0_2_V_DM1.npy')
command = wao.supervisor.rtc.get_command(0)
wao.supervisor.rtc.set_command(0,command)
command[88:] = -M2V_DM1@M_DM0_2_M_DM1[:,0