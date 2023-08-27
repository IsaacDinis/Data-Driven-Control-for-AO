#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/saxo_DM/compass/compass_param.py

from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
from matplotlib import pyplot as plt
from scipy.spatial import KDTree
import astropy.io.fits as pfits

pupil_diam = wao.supervisor.config.p_geom.get_pupdiam()
pupil_grid = make_pupil_grid(pupil_diam)
zernike_basis = make_zernike_basis(3, 1, pupil_grid)
tilt = zernike_basis[2].shaped
pupil_valid = zernike_basis[0].shaped

M2V_DM1 = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/M2V_DM1.npy')
M2V_DM0 = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/M2V_DM0.npy')
S2M_DM0 = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/S2M_DM0.npy')
S2M_DM1 = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/S2M_DM1.npy')
command = wao.supervisor.rtc.get_command(0)

M2V = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/M2V.npy')
M_DM0_2_M_DM1 = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/M_DM0_2_M_DM1.npy')

V2V = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/V_DM0_2_V_DM1.npy')
V2V2 = np.load('../../Data-Driven-Control-for-AO/saxo_DM/compass/calib_mat/V_DM1_2_V_DM0.npy')
# phase_tilt = pfits.getdata('../../Data-Driven-Control-for-AO/2DM_study/data3/dist_tilt.fits')

n_act_DM0 = wao.supervisor.config.p_dms[0].get_ntotact()
n_act_DM1 = wao.supervisor.config.p_dms[1].get_ntotact()


# pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data2/slopes4.fits", slopes, overwrite = True)
command = wao.supervisor.rtc.get_command(0)
wao.supervisor.rtc.set_command(0,command)
command[88:] = -M2V_DM1@M_DM0_2_M_DM1[:,0]*1000


mode_n = 1
amp = 1

u_DM0 = M2V_DM0[:,mode_n]
u_DM1 = M2V[:,n_act_DM0+mode_n]
# u_DM1 = M2V_DM1 @ M_DM0_2_M_DM1[:,mode_n]

command *=0 
command[:n_act_DM0] = u_DM0*amp*0
command[n_act_DM0:] = u_DM1[n_act_DM0:]*amp

wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()



mode_n = 1
amp = 1

u_DM0 = M2V_DM0[:,mode_n]
# u_DM1 = M2V_DM1@M_DM0_2_M_DM1[:,mode_n]
u_DM1 = V2V@u_DM0


command *=0 
command[:n_act_DM0] = u_DM0*amp
command[n_act_DM0:] = -u_DM1*amp

wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()


slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
print(modes_DM0[1])
print(modes_DM1[0])

mode_n = n_act_DM0
wao.supervisor.rtc.set_command(0, M2V[:,mode_n]) 
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()


# print(np.std(a))

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase,where = pupil_valid.astype(bool)))

print(np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid))

slopes = wao.supervisor.rtc.get_slopes(0)

modes_DM0 = np.dot(S2M_DM0,slopes)/557.2036425356519
modes_DM1 = np.dot(S2M_DM1,slopes)/557.2036425356519
print(modes_DM0[1])
print(modes_DM1[1])




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
command_HODM = -V2V@command[:88]
command_HODM[command_HODM>-0.0001] = 0
command_HODM[np.argmin(command_HODM)]=0
# command[88:] = command_HODM


wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
print(np.max(np.abs(command[88:])))
phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.max(np.abs(phase)))


HODM_act = 700
pos_LODM = np.array([wao.supervisor.config.p_dms[0].get_xpos(),wao.supervisor.config.p_dms[0].get_ypos()]).T
pos_HODM = np.array([wao.supervisor.config.p_dms[1].get_xpos(),wao.supervisor.config.p_dms[1].get_ypos()]).T
kd_tree_LODM = KDTree(pos_LODM)
kd_tree_HODM = KDTree(pos_HODM)

d, i = kd_tree_LODM.query(pos_HODM[HODM_act,:], k=4)
# d, i = kd_tree_LODM.query(np.array([1000,1000]), k=4)
w = 1/d
w /= np.sum(w)
# print(w, i, sep='\n')

command *= 0
command_LODM = np.zeros(88)
for act in range(4):
	command_LODM[i[act]] = w[act]/0.49920276
	command_HODM = -V2V@command_LODM
	command_HODM[command_HODM>-0.0001] = 0
	command += np.concatenate([command_LODM,command_HODM])
	command_LODM *= 0
command[88+HODM_act] = 0
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()

phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.max(np.abs(phase)))


command *= 0
command_LODM = np.zeros(88)
for act in range(4):
	command_LODM[i[act]] = w[act]#/0.43912143
command_HODM = -V2V@command_LODM
command_HODM[command_HODM>-0.0001] = 0
command += np.concatenate([command_LODM,command_HODM])
command[88+HODM_act] = 0
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()

phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.max(np.abs(phase)))



d, i = kd_tree_HODM.query(np.array([1000,1000]), k=4)
w = 1/d
w /= np.sum(w)
command_HODM[i] = w
command[88:] += command_HODM
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()


command[88+HODM_act] = -1/10
# command[:88] = -V2V2@command[88:]
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
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




m_pupil_size = wao.supervisor.get_m_pupil().shape[0]
pupil_grid = make_pupil_grid(m_pupil_size)
zernike_tel = make_zernike_basis(3, 1, pupil_grid)
tilt_tel = zernike_tel[1].shaped
pup_valid_tel = zernike_tel[0].shaped

# tilt_record = np.dstack([tilt]*phase_tilt.size)
tilt_record = np.dstack([tilt_tel]*10)
wao.supervisor.tel.set_input_phase(tilt_record)


slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
print(modes_DM0[1])
print(modes_DM1[1])


print(np.std(tilt,where = pupil_valid.astype(bool)))
print(np.max(np.abs(tilt*pupil_valid)))



wao.supervisor.tel.set_input_phase(tilt_record)
wao.supervisor.rtc.set_command(0,command*0)
wao.supervisor.next()
wao.supervisor.next()

# print(np.std(tilt,where = pupil_valid.astype(bool)))
tilt = zernike_basis[1].shaped
tip = zernike_basis[2].shaped

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase,where = pupil_valid.astype(bool)))
print(np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid))
print(np.sum(np.multiply(target_phase,tip))/np.sum(pupil_valid))
print(np.max(np.abs(target_phase*pupil_valid)))


print(np.std(tilt_record[:,:,0],where = pup_valid_tel.astype(bool)))
print(np.max(np.abs(tilt_record[:,:,0]*pup_valid_tel)))

mode_n = 1
# wao.supervisor.tel.reset_input_phase()
u_DM0 = M2V_DM0[:,mode_n]
u_DM1 = M2V_DM1[:,mode_n]
# u_DM1 = M2V_DM1 @ M_DM0_2_M_DM1[:,mode_n]

slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
print(modes_DM0[1])
print(modes_DM1[1])

# command *= 0
# command[:88] += -u_DM0*modes_DM0[mode_n]
# command[:88] -= M2V_DM0 @ modes_DM0
# command[88:] -= u_DM1*modes_DM0[1]
command[88:] -= M2V_DM1 @ modes_DM1

wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase,where = pupil_valid.astype(bool)))
print(np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid))
print(np.sum(np.multiply(target_phase,tip))/np.sum(pupil_valid))
print(np.max(np.abs(target_phase*pupil_valid)))

slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)*0.001763353
modes_DM1 = np.dot(S2M_DM1,slopes)*0.001763353

print(modes_DM0[mode_n])
print(modes_DM1[mode_n])

print(modes_DM0[2])
print(modes_DM1[2])









pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data3/DM0_closed_all.fits", modes_DM0, overwrite = True)
pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data3/DM1_closed_all.fits", modes_DM1, overwrite = True)





P2M_DM0 = pfits.getdata('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/P2M_DM0.fits')
P2M_DM1 = pfits.getdata('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/P2M_DM1.fits')
M2V_DM0 = pfits.getdata('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM0.fits')
M2V_DM1 = pfits.getdata('../../Data-Driven-Control-for-AO/2DM_study/compass/calib_mat/M2V_DM1.fits')
command = wao.supervisor.rtc.get_command(0)

mode_n = 1

u_DM0 = M2V_DM0[:,mode_n]
u_DM1 = M2V_DM1[:,mode_n]

# command *= 0
# command[:88] += -u_DM0*modes_DM0[1]
command[:88] -= M2V_DM0 @ modes_DM0
# command[88:] -= u_DM1*0.01
# command[88:] -= M2V_DM1 @ modes_DM1


wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()

phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
phase = phase[pupil_valid == 1]
modes_DM0 = np.dot(P2M_DM0,phase)
modes_DM1 = np.dot(P2M_DM1,phase)
print(modes_DM0[1])
print(modes_DM1[1])



target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase,where = pupil_valid.astype(bool)))
print(np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid))
print(np.sum(np.multiply(target_phase,tip))/np.sum(pupil_valid))
print(np.max(np.abs(target_phase*pupil_valid)))


pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data5/DM0_close_all.fits", modes_DM0, overwrite = True)
pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data5/DM1_close_all.fits", modes_DM1, overwrite = True)
pfits.writeto("../../Data-Driven-Control-for-AO/2DM_study/data5/phase.fits", phase, overwrite = True)




command[1] = (-1+0.02154211927797799)/9.546481186594894

wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()


target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase,where = pupil_valid.astype(bool)))
print(np.sum(np.multiply(target_phase,tilt))/np.sum(pupil_valid))
print(np.sum(np.multiply(target_phase,tip))/np.sum(pupil_valid))
print(np.max(np.abs(target_phase*pupil_valid)))