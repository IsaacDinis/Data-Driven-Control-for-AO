#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/2DM_study/compass/compass_param.py

from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 

pupil_diam = wao.supervisor.config.p_geom.get_pupdiam()
pupil_grid = make_pupil_grid(pupil_diam)
zernike_basis = make_zernike_basis(3, 1, pupil_grid)
tilt = zernike_basis[2].shaped
pupil_valid = zernike_basis[0].shaped

M2V_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM1.npy')
M2V_DM0 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM0.npy')
S2M_DM0 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/S2M_DM0.npy')
S2M_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/S2M_DM1.npy')
command = wao.supervisor.rtc.get_command(0)

M2V = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V.npy')
M_DM0_2_M_DM1 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M_DM0_2_M_DM1.npy')
nact_DM0 = 88
V2V = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM0_2_V_DM1.npy')
V2V2 = np.load('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/V_DM1_2_V_DM0.npy')

pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data2/slopes4.fits", slopes, overwrite = True)
command = wao.supervisor.rtc.get_command(0)
wao.supervisor.rtc.set_command(0,command)
command[88:] = -M2V_DM1@M_DM0_2_M_DM1[:,0]*1000


mode_n = 0
amp = 1.0

u_DM0 = M2V_DM0[:,mode_n]
u_DM1 = M2V_DM1[:,mode_n]
# u_DM1 = M2V_DM1 @ M_DM0_2_M_DM1[:,mode_n]

command[:88] = -u_DM0*amp
command[88:] = u_DM1*amp

wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()


# print(np.std(a))

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase,where = pupil_valid.astype(bool)))

print(np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid))

slopes = wao.supervisor.rtc.get_slopes(0)

modes_DM0 = np.dot(S2M_DM0,slopes)/557.2036425356519
modes_DM1 = np.dot(S2M_DM1,slopes)/557.2036425356519
print(modes_DM0[0])
print(modes_DM1[0])

u_DM0 = M2V_DM0[:,mode_n]
u_DM1 = M2V_DM1[:,mode_n]
# u_DM1 = M2V_DM1 @ M_DM0_2_M_DM1[:,mode_n] 

command[88:] = u_DM1*0
command[:88] = u_DM0*amp

wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()


# print(np.std(a))








target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(np.dot(target_phase,tilt),where = pupil_valid.astype(bool)))


slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)/2.1816739410448625
modes_DM1 = np.dot(S2M_DM1,slopes)/2.08673169704465
print(modes_DM0[0])
print(modes_DM1[0])


#DM1 1.6243452
#DM0 1.7792003

command[88:] = -u_DM1*amp/1.5526716
command[:88] = u_DM0*amp/1.5617224

0.034009222

-M2V_DM1 @ M_DM0_2_M_DM1 


u_DM0 = np.random.rand(88)*100
u_DM0 -= np.mean(u_DM0)

command[:88] = u_DM0
u_DM1 = -V2V@command[:88]
command[88:] = u_DM1-np.mean(u_DM1)


wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()

a = wao.supervisor.target.get_tar_phase(0,pupil=True)

print(np.std(a))



command *= 0

command[500] = 0.5 
# command[:88] = -V2V2@command[88:]
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()


command *= 0

command[40] = 0.1 
command[88:] = -V2V@command[:88]

command[88:] -= np.mean(command[88:])
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
print(np.max(np.abs(command[88:])))
phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.max(np.abs(phase)))


command *= 0
command[500] = 0.5 
command[499] = 0.5 
command[460] = 0.5 
command[461] = 0.5 
command[25] = -0.75
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.max(np.abs(phase)))