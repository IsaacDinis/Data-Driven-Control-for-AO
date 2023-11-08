#ipython -i shesha/widgets/widget_ao.py ~/Data-Driven-Control-for-AO/ristretto_bump/compass/ristretto_param.py

from hcipy.field import make_pupil_grid 
from hcipy.mode_basis import make_zernike_basis 
from matplotlib import pyplot as plt
from scipy.spatial import KDTree
import astropy.io.fits as pfits
from scipy.spatial import KDTree
from scipy.interpolate import interpn
import utils
# print(wao.supervisor.config.p_dms[0].get_ntotact())
# print(wao.supervisor.config.p_dms[1].get_ntotact())
# print(wao.supervisor.rtc.get_slopes(0).shape)

M2V_DM1 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/M2V_DM1.npy')
M2V_DM0 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/M2V_DM0.npy')
S2M_DM0 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/S2M_DM0.npy')
S2M_DM1 = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/S2M_DM1.npy')
command = wao.supervisor.rtc.get_command(0)
M2V = np.load('../../Data-Driven-Control-for-AO/ristretto/compass/calib_mat/M2V.npy')

n_act0 = wao.supervisor.config.p_dms[0].get_ntotact()
n_act1 = wao.supervisor.config.p_dms[1].get_ntotact()
n_act2 = wao.supervisor.config.p_dms[2].get_ntotact()

n_modes = 500

wao.supervisor.rtc.set_command(0,np.concatenate((M2V_DM0[:,0]*0, M2V_DM1[:,n_modes]), axis=0))
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

wao.supervisor.rtc.set_command(0,M2V[:,n_act0+800]*0.01)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()
slopes = wao.supervisor.rtc.get_slopes(0)

modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
print(modes_DM0[1])
print(modes_DM1[1])

wao.supervisor.rtc.set_command(0,M2V[:,0]*0)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

wao.supervisor.rtc.set_command(0,np.concatenate((M2V_DM0[:,0]*0, M2V_DM1[:,n_modes]*0), axis=0))
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

n_modes = 500
slopes = wao.supervisor.rtc.get_slopes(0)
modes_DM0 = np.dot(S2M_DM0,slopes)
modes_DM1 = np.dot(S2M_DM1,slopes)
wao.supervisor.rtc.set_command(0,np.concatenate((M2V_DM0[:,0]*0, M2V_DM1[:,:n_modes]@modes_DM1[:n_modes]), axis=0))
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()



print(modes_DM0[1])
print(modes_DM1[1])

norm = np.zeros(800)
for i in range(800):
    wao.supervisor.rtc.set_command(0,M2V[:,n_act0+i])
    wao.supervisor.next()
    wao.supervisor.next()
    wao.supervisor.next()
    target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
    norm[i] = np.std(target_phase)
norm /=np.max(norm)   
plt.plot(norm)
plt.show()



pos_LODM = np.array([wao.supervisor.config.p_dms[0].get_xpos(),wao.supervisor.config.p_dms[0].get_ypos()]).T
pos_HODM = np.array([wao.supervisor.config.p_dms[1].get_xpos(),wao.supervisor.config.p_dms[1].get_ypos()]).T
pos_HODM -= np.min(pos_HODM,axis = 0)

command *= 0
command[-1] = -1
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
plt.imshow(target_phase)


# 685 686 606 765
p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

middle = (np.max(xpos)-np.min(xpos))/2
# xpos0 = 240.5 # 960 961
# ypos0 = 184.5
xpos0 = 240.0 # 942 942 980 981
ypos0 = 185.8

plt.scatter(xpos, ypos, marker='.', color="red")
plt.scatter(xpos0, ypos0, marker='.', color="green")

influ = p_dm._influ
plt.imshow(influ[:,:,800])
influ0 = p_dm._influ[:,:,0]
influ0 = np.expand_dims(influ0, axis=2)
i10 = xpos0-20.5
j10 = ypos0-20.5 
xcenter = p_geom.cent
ycenter = p_geom.cent

xpos0 += p_dm._n1
ypos0 += p_dm._n1
i10 += p_dm._n1
j10 += p_dm._n1
file_name = 'bump.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i10,j10,influ0,xpos0,ypos0,xcenter,ycenter,pixsize,diam)





# 685 686 606 765
p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

i1 = p_dm._i1.copy()
j1 = p_dm._j1.copy()
influ = p_dm._influ.copy()

middle = (np.max(xpos)-np.min(xpos))/2
# xpos0 = 240.5 # 960 961
# ypos0 = 184.5
xpos0 = 240.0 # 960 961
ypos0 = 185.8

plt.scatter(xpos, ypos, marker='.', color="red")
plt.scatter(xpos0, ypos0, marker='.', color="green")
plt.scatter(xpos[960], ypos[960], marker='.', color="blue")


plt.imshow(influ[:,:,800])
influ0 = p_dm._influ[:,:,0]
influ0 = np.expand_dims(influ0, axis=2)

i10 = xpos0-20.5
j10 = ypos0-20.5 

p_dm._i1
xcenter = p_geom.cent
ycenter = p_geom.cent

xpos = np.append(xpos,xpos0)
ypos = np.append(ypos,ypos0)

i1 = np.append(i1,i10)
j1 = np.append(j1,j10)

influ = np.append(influ,influ0,axis = 2)

xpos += p_dm._n1
ypos += p_dm._n1
i1 += p_dm._n1
j1 += p_dm._n1

file_name = 'bump.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i1,j1,influ,xpos,ypos,xcenter,ycenter,pixsize,diam)





p_dm._i1
wao.supervisor.config.p_dms[2].set_n1(331)


HODM_act = 700
pos_LODM = np.array([wao.supervisor.config.p_dms[0].get_xpos(),wao.supervisor.config.p_dms[0].get_ypos()]).T
pos_HODM = np.array([wao.supervisor.config.p_dms[2].get_xpos(),wao.supervisor.config.p_dms[2].get_ypos()]).T
kd_tree_LODM = KDTree(pos_LODM)


d, i = kd_tree_LODM.query(pos_HODM[950,:], k=4)
w = 1/d
w /= np.sum(w)

command *= 0
# command[-1] = -1
command[n_act0+n_act1+950] = -1
command_LODM = np.zeros(n_act0)
# command_dead_act = np.zeros(n_actus_DM0 + n_actus_DM1)
for act in range(4):
    command_LODM[i[act]] = w[act]/0.61707693
    # command_HODM = -V_DM0_2_V_DM1@command_LODM
    # command_HODM[command_HODM>-0.0001] = 0
    # command_dead_act += np.concatenate([command_LODM,command_HODM])
    command[:n_act0] += command_LODM
    command_LODM *= 0
    # command_dead_act[n_actus_DM0+HODM_act] = 0


wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
print(np.std(target_phase))
print(np.max(target_phase))


pos_HODM = np.array([wao.supervisor.config.p_dms[1].get_xpos(),wao.supervisor.config.p_dms[1].get_ypos()]).T
p_dm =  wao.supervisor.config.p_dms[1]
p_geom = wao.supervisor.config.p_geom

command *= 0
command[n_act0+941] = 3.5
command[n_act0+942] = 3.5
command[n_act0+980] = 3.5
command[n_act0+981] = 3.5
command[-1] = -np.mean([command[n_act0+941],command[n_act0+942],command[n_act0+980],command[n_act0+981]])*1.7
# command[-1] = -7
wao.supervisor.rtc.set_command(0,command)
wao.supervisor.next()
wao.supervisor.next()
wao.supervisor.next()
target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)




plt.imshow(target_phase)
plt.scatter(pos_HODM[942,1]-p_geom._p1, pos_HODM[942,0]-p_geom._p1, marker='.', color="blue")
plt.scatter(pos_HODM[941,1]-p_geom._p1, pos_HODM[941,0]-p_geom._p1, marker='.', color="blue")
plt.scatter(pos_HODM[980,1]-p_geom._p1, pos_HODM[980,0]-p_geom._p1, marker='.', color="blue")
plt.scatter(pos_HODM[981,1]-p_geom._p1, pos_HODM[981,0]-p_geom._p1, marker='.', color="blue")

m_l = 145
m_r = 185
m_u = 200
m_d = 240

pos_act = np.array([[pos_HODM[942,1]-p_geom._p1-m_l,pos_HODM[942,0]-p_geom._p1-m_u],
    [pos_HODM[941,1]-p_geom._p1-m_l,pos_HODM[941,0]-p_geom._p1-m_u],
    [pos_HODM[980,1]-p_geom._p1-m_l,pos_HODM[980,0]-p_geom._p1-m_u],
    [pos_HODM[981,1]-p_geom._p1-m_l,pos_HODM[981,0]-p_geom._p1-m_u]])

pfits.writeto('../..//Data-Driven-Control-for-AO/ristretto_bump/matlab/pos_act.fits', pos_act, overwrite = True)
target_phase_croped = target_phase[m_u:m_d,m_l:m_r]
pfits.writeto('../..//Data-Driven-Control-for-AO/ristretto_bump/matlab/target_phase_croped_35.fits', target_phase_croped, overwrite = True)

plt.imshow(target_phase_croped)
plt.scatter(pos_HODM[942,1]-p_geom._p1-m_l, pos_HODM[942,0]-p_geom._p1-m_u, marker='.', color="blue")
plt.scatter(pos_HODM[941,1]-p_geom._p1-m_l, pos_HODM[941,0]-p_geom._p1-m_u, marker='.', color="blue")
plt.scatter(pos_HODM[980,1]-p_geom._p1-m_l, pos_HODM[980,0]-p_geom._p1-m_u, marker='.', color="blue")
plt.scatter(pos_HODM[981,1]-p_geom._p1-m_l, pos_HODM[981,0]-p_geom._p1-m_u, marker='.', color="blue")

m_size = wao.supervisor.config.p_geom.get_mpupil().shape[0]
s_size = wao.supervisor.config.p_geom.get_spupil().shape[0]
pad_size = int((m_size-s_size)/2)

flat *= 1e6

flat_size = flat.shape[0]
dummy_x,dummy_y = np.meshgrid(np.linspace(0,flat_size-10,s_size),np.linspace(0,flat_size-10,s_size))


flat_interp = interpn((np.arange(flat_size), np.arange(flat_size)), flat,(dummy_x, dummy_y ))

flat_pad = np.pad(flat_interp,pad_size)
flat_record = np.dstack([flat_pad])
wao.supervisor.tel.set_input_phase(flat_record)

plt.imshow(flat)

plt.imshow(flat_pad)


flat_record = np.dstack([flat])
wao.supervisor.tel.set_input_phase(flat_record)


flat = pfits.getdata('../../Data-Driven-Control-for-AO/ristretto_true_DM/compass/calib_mat/flat.fits')
wao.supervisor.tel.set_input_phase(flat)
plt.imshow(flat)

target_phase = wao.supervisor.target.get_tar_phase(0,pupil=True)
plt.imshow(target_phase)



act_pos_43 = pfits.getdata('../../Data-Driven-Control-for-AO/ristretto_true_DM/compass/calib_mat/act_pos.fits')

p_geom = wao.supervisor.config.p_geom
p_tel =  wao.supervisor.config.p_tel
p_dm =  wao.supervisor.config.p_dms[1]
pixsize = wao.supervisor.config.p_geom.get_pixsize()
diam = p_tel.get_diam()
xpos = p_dm._xpos-p_dm._n1
ypos = p_dm._ypos-p_dm._n1

i1 = p_dm._i1.copy()
j1 = p_dm._j1.copy()
influ = p_dm._influ.copy()



min_xpos = np.min(xpos)
min_ypos = np.min(ypos)

xpos43 = act_pos_43[:,0]-np.min(act_pos_43[:,0])
ypos43 = act_pos_43[:,1]-np.min(act_pos_43[:,1])
xpos -= min_xpos
ypos -= min_ypos

xpos43 /= np.max(xpos43)/np.max(xpos)
ypos43 /= np.max(ypos43)/np.max(ypos)

xpos43 += min_xpos
ypos43 += min_ypos
xpos += min_xpos
ypos += min_ypos

plt.scatter(xpos, ypos, marker='.', color="red")
# plt.scatter(xpos0, ypos0, marker='.', color="green")
plt.scatter(xpos43, ypos43, marker='.', color="blue")


plt.imshow(influ[:,:,800])


i1 = xpos43 - 16.5
j1 = ypos43 - 16.5


xcenter = p_geom.cent
ycenter = p_geom.cent





influ = np.append(influ,influ[:,:,:xpos43.shape[0]-xpos.shape[0]],axis = 2)

xpos43 += p_dm._n1
ypos43 += p_dm._n1
xpos += p_dm._n1
ypos += p_dm._n1
i1 += p_dm._n1
j1 += p_dm._n1

file_name = 'ristretto_43.fits'
dm_custom = utils.write_dm_custom_fits(file_name,i1,j1,influ,xpos43,ypos43,xcenter,ycenter,pixsize,diam)


