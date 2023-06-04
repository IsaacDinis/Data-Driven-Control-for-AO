M2V_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM1.npy')
M2V_DM0 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM0.npy')
S2M_DM0 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/S2M_DM0.npy')
S2M_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/S2M_DM1.npy')
command = wao.supervisor.rtc.get_command(0)

M2V = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V.npy')
M_DM0_2_M_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M_DM0_2_M_DM1.npy')
nact_DM0 = 88
V2V = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM0_2_V_DM1.npy')
pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data2/slopes4.fits", slopes, overwrite = True)
command = wao.supervisor.rtc.get_command(0)
wao.supervisor.rtc.set_command(0,command)
command[88:] = -M2V_DM1@M_DM0_2_M_DM1[:,0]


mode_n = 0
amp = 1

u_DM0 = M2V_DM0[:,mode_n]
u_DM1 = M2V_DM1[:,mode_n]

command[88:] = u_DM1*amp
command[:88] = u_DM0*0
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()

a = wao.supervisor.target.get_tar_phase(0)*1000

print(np.std(a))




slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
print(modes_DM0[0])
print(modes_DM1[0])